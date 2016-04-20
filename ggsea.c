
/**
 * @file ggsea.c
 *
 * @brief Graph-to-Graph Seed-and-Extend Alignment
 *
 * @author Hajime Suzuki
 * @date 2016/4/12
 */

#include <stdint.h>
#include "ggsea.h"
#include "hmap/hmap.h"
#include "gref/gref.h"
#include "gaba/gaba.h"
#include "psort/psort.h"
#include "kvec.h"
#include "sassert.h"
#include "log.h"

#define UNITTEST_UNIQUE_ID			10
#define UNITTEST 					1

#include "unittest.h"

/* inline directive */
#define _force_inline				inline

/* constants */
#define MARGIN_SEQ_SIZE				( 64 )
#define MARGIN_SEQ_OFFSET			( 16 )
#define MARGIN_SEQ_LEN				( 32 )

/* assertions for ggsea_params_s */
_static_assert_offset(struct ggsea_params_s, seq_a_format, struct gaba_params_s, seq_a_format, 0);
_static_assert_offset(struct ggsea_params_s, seq_a_direction, struct gaba_params_s, seq_a_direction, 0);
_static_assert_offset(struct ggsea_params_s, seq_b_format, struct gaba_params_s, seq_b_format, 0);
_static_assert_offset(struct ggsea_params_s, seq_b_direction, struct gaba_params_s, seq_b_direction, 0);
_static_assert_offset(struct ggsea_params_s, xdrop, struct gaba_params_s, xdrop, 0);
_static_assert_offset(struct ggsea_params_s, score_matrix, struct gaba_params_s, score_matrix, 0);

/**
 * @struct ggsea_conf_s
 * @brief configuration container
 */
struct ggsea_conf_s {
	/* alignment context */
	gaba_t *gaba;

	/* params */
	int64_t rep_kmer_hmap_size;
	int64_t max_rep_vec_size;		/* max kmer vector size */
	struct ggsea_params_s params;
};

/**
 * @fn ggsea_conf_init
 */
ggsea_conf_t *ggsea_conf_init(
	ggsea_params_t const *params)
{
	/* restore default params */
	struct ggsea_params_s const default_params = { 0 };
	struct ggsea_params_s p = (params == NULL) ? default_params : *params;

	/* restore defaults */
	#define restore(param, def)		{ (param) = ((uint64_t)(param) == 0) ? (def) : (param); }

	restore(p.kmer_cnt_thresh, 100);
	restore(p.popcnt_thresh, 25);
	restore(p.popcnt_low_thresh, 15);
	restore(p.overlap_thresh, 3);
	restore(p.num_threads, 0);

	#undef restore

	/* malloc mem */
	struct ggsea_conf_s *conf = (struct ggsea_conf_s *)malloc(
		sizeof(struct ggsea_conf_s));
	if(conf == NULL) {
		return(NULL);
	}
	memset(conf, 0, sizeof(struct ggsea_conf_s));

	/* init alignment context */
	conf->gaba = gaba_init((struct gaba_params_s const *)&p);
	if(conf->gaba == NULL) {
		free(conf);
		return(NULL);
	}

	/* store constants */
	conf->rep_kmer_hmap_size = 1024;
	conf->max_rep_vec_size = 128;
	conf->params = p;

	return((ggsea_conf_t *)conf);
}

/**
 * @fn ggsea_conf_clean
 */
void ggsea_conf_clean(
	ggsea_conf_t *_conf)
{
	struct ggsea_conf_s *conf = (struct ggsea_conf_s *)_conf;

	if(conf != NULL) {
		gaba_clean(conf->gaba); conf->gaba = NULL;
		free(conf);
	}
	return;
}

/**
 * @struct ggsea_rep_kmer_cont_s
 */
struct ggsea_rep_kmer_cont_s {
	hmap_header_t header;
	kvec_t(struct gref_gid_pos_s) rv;
	kvec_t(struct gref_gid_pos_s) qv;
};

/**
 * @struct ggsea_weak_kmer_s
 */
struct ggsea_weak_kmer_s {
	// uint64_t kmer;			/* kmer info is not necessary */
	struct gref_gid_pos_s rpos;
	struct gref_gid_pos_s qpos;
};

/**
 * @struct ggsea_node_s
 * @brief element of node array
 */
struct ggsea_node_s {
	uint64_t incoming_edges;
};

/**
 * @struct ggsea_segq_s
 * @brief segment info container (to push into heapqueue)
 */
struct ggsea_segq_s {
	gaba_fill_t const *fill;
	uint32_t rgid;
	uint32_t qgid;
};

/**
 * @struct ggsea_fill_pair_s
 */
struct ggsea_fill_pair_s {
	gaba_fill_t const *fw;
	gaba_fill_t const *rv;
};

/**
 * @struct ggsea_ctx_s
 */
struct ggsea_ctx_s {
	/* constants */
	struct ggsea_conf_s conf;	/* copy of conf */

	/* current seqeunce info */
	gref_idx_t const *r;
	gref_acv_t const *q;

	/* repetitive kmer container */
	hmap_t *rep_kmer;

	/* weak kmer list */
	kvec_t(struct ggsea_weak_kmer_s) weak_kmer;

	/* dp context */
	gaba_dp_t *dp;
	struct ggsea_node_s *node;			/* node info array */
	kvec_t(struct ggsea_segq_s) dagq;	/* queue with segments in dag subgraph */
	uint8_t *margin;
	struct gaba_section_s fw_margin, rv_margin;

	/* max */
	gaba_fill_t const *max, *fw_max, *rv_max;

	/* score order info */

	/* red-black-tree に入れるんだったっけ?? */

	/* result container */
	kvec_t(struct ggsea_fill_pair_s) res;
};

/**
 * @fn ggsea_ctx_init
 */
ggsea_ctx_t *ggsea_ctx_init(
	ggsea_conf_t const *conf,
	gref_idx_t const *ref)
{
	/* malloc mem */
	struct ggsea_ctx_s *ctx = (struct ggsea_ctx_s *)malloc(
		sizeof(struct ggsea_ctx_s));
	if(ctx == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	memset(ctx, 0, sizeof(struct ggsea_ctx_s));

	/* init contexts */
	ctx->rep_kmer = hmap_init(
		conf->rep_kmer_hmap_size, sizeof(struct ggsea_rep_kmer_cont_s));
	if(ctx->rep_kmer == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	/* copy conf */
	ctx->conf = *conf;

	/* clear ref and query pointers */
	if(ref == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	ctx->r = ref;
	ctx->q = NULL;

	/* init dp context (with invalid seq pair) */
	ctx->dp = gaba_dp_init(conf->gaba,
		&gaba_build_seq_pair(gref_get_ptr(ref), gref_get_total_len(ref), NULL, 0));
	if(ctx->dp == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	/* init node array */
	uint32_t sec_cnt = gref_get_section_cnt(ref);
	ctx->node = (struct ggsea_node_s *)malloc(sizeof(struct ggsea_node_s) * sec_cnt);
	if(ctx->node == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	/* init queue */
	kv_hq_init(ctx->dagq);
	if(kv_ptr(ctx->dagq) == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	/* init margin seq */
	ctx->margin = (uint8_t *)malloc(2 * sizeof(uint8_t) * MARGIN_SEQ_SIZE);
	if(ctx->margin == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	uint8_t const pad[7] = { [GABA_ASCII] = 'A' };

	/* init forward margin section */
	memset(ctx->margin,
		pad[ctx->conf.params.seq_a_format],
		MARGIN_SEQ_SIZE);
	ctx->fw_margin = (struct gaba_section_s){
		.id = 0xfffc,
		.len = MARGIN_SEQ_LEN,
		.base = ctx->margin - gref_get_ptr(ctx->r)
	};

	/* init reverse margin section */
	memset(&ctx->margin[MARGIN_SEQ_SIZE],
		pad[ctx->conf.params.seq_b_format],
		MARGIN_SEQ_SIZE);
	ctx->rv_margin = (struct gaba_section_s){
		.id = 0xfffd,
		.len = MARGIN_SEQ_LEN,
		.base = &ctx->margin[MARGIN_SEQ_SIZE] - gref_get_ptr(ctx->r)
	};

	/* init max */
	ctx->max = ctx->fw_max = ctx->rv_max = NULL;

	/* init result array */
	kv_init(ctx->res);

	return(ctx);

_ggsea_ctx_init_error_handler:;
	if(ctx != NULL) {
		hmap_clean(ctx->rep_kmer); ctx->rep_kmer = NULL;
		gaba_dp_clean(ctx->dp); ctx->dp = NULL;
		free(ctx->node); ctx->node = NULL;
		kv_hq_destroy(ctx->dagq);
		free(ctx->margin); ctx->margin = NULL;
		kv_destroy(ctx->res);
		free(ctx);
	}
	return(NULL);
}

/**
 * @fn ggsea_ctx_clean
 */
void ggsea_ctx_clean(
	ggsea_ctx_t *ctx)
{
	if(ctx != NULL) {
		hmap_clean(ctx->rep_kmer); ctx->rep_kmer = NULL;
		gaba_dp_clean(ctx->dp); ctx->dp = NULL;
		free(ctx->node); ctx->node = NULL;
		kv_hq_destroy(ctx->dagq);
		free(ctx->margin); ctx->margin = NULL;
		kv_destroy(ctx->res);
		free(ctx);
	}
	return;
}

/**
 * @fn ggsea_ctx_flush
 */
static _force_inline
void ggsea_ctx_flush(
	struct ggsea_ctx_s *ctx,
	gref_acv_t const *query)
{
	/* set sequence info */
	ctx->q = query;

	/* flush queues */
	kv_hq_clear(ctx->dagq);

	/* flush dp context for the new read */
	gaba_dp_flush(ctx->dp,
		&gaba_build_seq_pair(
			gref_get_ptr(ctx->r), gref_get_total_len(ctx->r),
			gref_get_ptr(ctx->q), gref_get_total_len(ctx->q)));

	/* flush result array */
	kv_clear(ctx->res);
	return;
}

/**
 * @fn ggsea_dedup_rep_kmer
 */
static _force_inline
int64_t ggsea_dedup_rep_kmer(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s *ptr,
	int64_t size)
{
	/* sort */
	psort_full((void *)ptr, size, sizeof(struct gref_gid_pos_s), ctx->conf.params.num_threads);

	/* dedup */
	uint64_t *p = (uint64_t *)ptr;
	int64_t cnt = 0;
	for(int64_t i = 0; i < size; i++) {
		if(p[cnt] == p[i]) { continue; }
		p[++cnt] = p[i];
	}
	return(cnt);
}

/**
 * @fn ggsea_save_rep_kmer
 */
static _force_inline
void ggsea_save_rep_kmer(
	struct ggsea_ctx_s *ctx,
	uint64_t kmer,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* add kmer to the hashmap */
	uint32_t id = hmap_get_id(ctx->rep_kmer,
		(char const *)&kmer, sizeof(uint64_t));
	struct ggsea_rep_kmer_cont_s *c = hmap_get_object(ctx->rep_kmer, id);
	
	/* add rpos */
	kv_push(c->rv, rpos);
	if(kv_size(c->rv) > ctx->conf.max_rep_vec_size) {
		kv_size(c->rv) = ggsea_dedup_rep_kmer(ctx, kv_ptr(c->rv), kv_size(c->rv));
	}

	kv_push(c->qv, qpos);
	if(kv_size(c->qv) > ctx->conf.max_rep_vec_size) {
		kv_size(c->qv) = ggsea_dedup_rep_kmer(ctx, kv_ptr(c->qv), kv_size(c->qv));
	}
	return;
}

/**
 * @fn ggsea_save_weak_kmer
 */
static _force_inline
void ggsea_save_weak_kmer(
	struct ggsea_ctx_s *ctx,
	int64_t popcnt_score,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* discard seed if popcnt_score is below the threshold */
	if(popcnt_score < ctx->conf.params.popcnt_low_thresh) {
		return;
	}

	/* push to vector */
	kv_push(ctx->weak_kmer, ((struct ggsea_weak_kmer_s){
		.rpos = rpos,
		.qpos = qpos
	}));
	return;
}

/**
 * @fn ggsea_popcnt_filter
 */
static _force_inline
int64_t ggsea_popcnt_filter(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	return(0);
}

/**
 * @fn ggsea_update_link_mask
 */
static _force_inline
uint64_t ggsea_update_link_mask(
	struct ggsea_ctx_s *ctx,
	uint32_t rgid,
	uint32_t qgid)
{
	return(0);
}

/**
 * @fn ggsea_get_full_mask
 */
static _force_inline
uint64_t ggsea_get_full_mask(
	struct ggsea_ctx_s *ctx,
	uint32_t gid)
{
	return(0);
}

/**
 * @fn ggsea_get_ndag_section
 */
static _force_inline
struct ggsea_segq_s ggsea_get_ndag_section(
	struct ggsea_ctx_s *ctx)
{
	return((struct ggsea_segq_s){ 0 });
}

/**
 * @fn ggsea_extend_update_r
 */
static _force_inline
void ggsea_extend_update_r(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	/* update rsec */
	struct gref_link_s rlink = gref_get_link(ctx->r, rsec->gid);

	/* if len == 0, push margin */
	if(rlink.len == 0) {
		;
	}

	/* push section pair if all the incoming edges are scaned */
	for(int64_t i = 0; i < rlink.len; i++) {
		uint32_t rgid = rlink.gid_arr[i];

		/* update link mask */
		uint64_t mask = ggsea_update_link_mask(ctx, rgid, rsec->gid);

		/* get full mask */
		uint64_t full = ggsea_get_full_mask(ctx, rgid);

		if(mask != full) { continue; }

		/* push queue */
		kv_hq_push(ctx->dagq, ((struct ggsea_segq_s){
			.fill = fill,
			.rgid = rgid,
			.qgid = qsec->gid
		}));
	}
	return;
}

/**
 * @fn ggsea_extend_update_q
 */
static _force_inline
void ggsea_extend_update_q(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	/* update qsec */
	struct gref_link_s qlink = gref_get_link(ctx->q, qsec->gid);

	/* if len == 0, push margin */
	if(qlink.len == 0) {
		;
	}

	/* push section pair if all the incoming edges are scaned */
	for(int64_t i = 0; i < qlink.len; i++) {
		uint32_t qgid = qlink.gid_arr[i];

		/* update link mask */
		uint64_t mask = ggsea_update_link_mask(ctx, qgid, qsec->gid);

		/* get full mask */
		uint64_t full = ggsea_get_full_mask(ctx, qgid);

		if(mask != full) { continue; }

		/* push queue */
		kv_hq_push(ctx->dagq, ((struct ggsea_segq_s){
			.fill = fill,
			.rgid = rsec->gid,
			.qgid = qgid
		}));
	}
	return;
}

/**
 * @fn ggsea_extend_update_rq
 */
static _force_inline
void ggsea_extend_update_rq(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	/* update sections */
	struct gref_link_s rlink = gref_get_link(ctx->r, rsec->gid);
	struct gref_link_s qlink = gref_get_link(ctx->q, qsec->gid);

	if((rlink.len | qlink.len) == 0) {
		;
	}

	/* push section pair if all the incoming edges are scaned */
	for(int64_t i = 0; i < rlink.len; i++) {
		uint32_t rgid = rlink.gid_arr[i];

		/* update link mask */
		uint64_t rmask = ggsea_update_link_mask(ctx, rgid, rsec->gid);

		/* get full mask */
		uint64_t rfull = ggsea_get_full_mask(ctx, rgid);

		if(rmask != rfull) { continue; }

		for(int64_t j = 0; j < qlink.len; j++) {
			uint32_t qgid = qlink.gid_arr[j];

			/* update link mask */
			uint64_t qmask = ggsea_update_link_mask(ctx, rgid, rsec->gid);

			/* get full mask */
			uint64_t qfull = ggsea_get_full_mask(ctx, rgid);

			if(qmask != qfull) { continue; }

			/* push queue */
			kv_hq_push(ctx->dagq, ((struct ggsea_segq_s){
				.fill = fill,
				.rgid = rgid,
				.qgid = qgid
			}));
		}
	}
	return;
}

/**
 * @fn ggsea_extend_update
 */
static _force_inline
void ggsea_extend_update(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	uint64_t rup = (fill->status & GABA_STATUS_UPDATE_A) != 0;
	uint64_t qup = (fill->status & GABA_STATUS_UPDATE_B) != 0;

	if(rup && qup) {
		ggsea_extend_update_rq(ctx, fill, rsec, qsec);
	} else if(rup) {
		ggsea_extend_update_r(ctx, fill, rsec, qsec);
	} else if(qup) {
		ggsea_extend_update_q(ctx, fill, rsec, qsec);
	} else {
		/* never reaches here */
	}
	return;
}

/**
 * @fn ggsea_extend_loop
 */
static _force_inline
int ggsea_extend_loop(
	struct ggsea_ctx_s *ctx)
{
	/* pop next section */
	struct ggsea_segq_s seg = (kv_hq_size(ctx->dagq) > 0)
		? kv_hq_pop(ctx->dagq)
		: ggsea_get_ndag_section(ctx);

	if(seg.fill == NULL) {
		return(-1);
	}

	/* extend */
	struct gref_section_s const *rsec = gref_get_section(ctx->r, seg.rgid);
	struct gref_section_s const *qsec = gref_get_section(ctx->q, seg.qgid);
	gaba_fill_t *fill = gaba_dp_fill(ctx->dp,
		seg.fill,
		(struct gaba_section_s *)rsec,
		(struct gaba_section_s *)qsec);

	/* update max */
	ctx->max = (fill->max > ctx->fw_max->max) ? fill : ctx->max;

	/* check xdrop term */
	if(fill->status & GABA_STATUS_TERM) {
		return(0);
	}

	/* check update flag */
	ggsea_extend_update(ctx, fill, rsec, qsec);
	return(0);
}

/**
 * @fn ggsea_extend
 */
static _force_inline
void ggsea_extend(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* initialize forward root */
	struct gref_section_s const *rsec = gref_get_section(ctx->r, rpos.gid);
	struct gref_section_s const *qsec = gref_get_section(ctx->q, qpos.gid);
	gaba_fill_t const *fill = ctx->max = gaba_dp_fill_root(ctx->dp,
		(struct gaba_section_s *)rsec, rpos.pos,
		(struct gaba_section_s *)qsec, qpos.pos);

	/* check xdrop term */
	if((fill->status & GABA_STATUS_TERM) == 0) {
		/* update first joint */
		ggsea_extend_update(ctx, fill, rsec, qsec);

		/* loop */
		while(ggsea_extend_loop(ctx) == 0) {}
	}

	/* save max */
	ctx->fw_max = ctx->max;

	/* initialize reverse root */
	rsec = gref_get_section(ctx->r, gref_rev(rpos.gid));
	qsec = gref_get_section(ctx->q, gref_rev(qpos.gid));
	fill = ctx->max = gaba_dp_fill_root(ctx->dp,
		(struct gaba_section_s *)rsec, rpos.pos,
		(struct gaba_section_s *)qsec, qpos.pos);	

	/* chech xdrop term */
	if((fill->status & GABA_STATUS_TERM) == 0) {
		/* update first joint */
		ggsea_extend_update(ctx, fill, rsec, qsec);

		/* loop */
		while(ggsea_extend_loop(ctx) == 0) {}
	}

	/* save max */
	ctx->rv_max = ctx->max;

	/* return max pair */
	kv_push(ctx->res, ((struct ggsea_fill_pair_s){
		.fw = ctx->fw_max,
		.rv = ctx->rv_max
	}));
}

/**
 * @fn ggsea_refine_results
 */
static _force_inline
ggsea_result_t ggsea_refine_results(
	struct ggsea_ctx_s *ctx)
{
	/* traceback */
	kvec_t(struct gaba_result_s const *) aln;
	kv_reserve(aln, kv_size(ctx->res));
	for(int64_t i = 0; i < kv_size(ctx->res); i++) {
		struct gaba_result_s const *r = gaba_dp_trace(ctx->dp,
			kv_at(ctx->res, i).fw, kv_at(ctx->res, i).rv, NULL);
		kv_push(aln, r);
	}

	/* pack result */
	struct ggsea_result_s r = {
		.aln = kv_ptr(aln),
		.cnt = kv_size(aln)
	};

	/* finish */
	return((ggsea_result_t)r);
}

/**
 * @fn ggsea_overlap_filter
 */
static _force_inline
int64_t ggsea_overlap_filter(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	int64_t score = 0;
	return(score);
}

/**
 * @fn ggsea_save_overlap_kmer
 */
static _force_inline
void ggsea_save_overlap_kmer(
	struct ggsea_ctx_s *ctx,
	int64_t score,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	return;
}

/**
 * @fn ggsea_finish
 */
static _force_inline
int ggsea_finish(
	struct ggsea_ctx_s *ctx)
{
	return(1);
}

/**
 * @fn ggsea_align
 */
ggsea_result_t ggsea_align(
	ggsea_ctx_t *_ctx,
	gref_acv_t const *query)
{
	struct ggsea_ctx_s *ctx = (struct ggsea_ctx_s *)_ctx;

	/* flush ctxing buffer */
	ggsea_ctx_flush(ctx, query);

	/* seeding */
	gref_iter_t *iter = gref_iter_init(ctx->q);

	struct gref_kmer_tuple_s t;
	while((t = gref_iter_next(iter)).kmer != GREF_ITER_KMER_TERM) {
		
		/* match */
		struct gref_match_res_s m = gref_match_2bitpacked(ctx->r, t.kmer);

		/* check if kmer is repetitive */
		if(m.len > ctx->conf.params.kmer_cnt_thresh) {
			ggsea_save_rep_kmer(ctx, t.kmer, t.gid_pos, m.gid_pos_arr[0]);
			continue;
		}

		/* for all positions on ref matched with the kmer */
		for(int64_t i = 0; i < m.len; i++) {

			/* first apply popcnt filter */
			int64_t popcnt_score = ggsea_popcnt_filter(ctx, t.gid_pos, m.gid_pos_arr[i]);
			if(popcnt_score <= ctx->conf.params.popcnt_thresh) {
				ggsea_save_weak_kmer(ctx, popcnt_score, t.gid_pos, m.gid_pos_arr[i]);
				continue;
			}

			/* passed, apply overlap filter */
			int64_t overlap_score = ggsea_overlap_filter(ctx, t.gid_pos, m.gid_pos_arr[i]);
			if(overlap_score <= ctx->conf.params.overlap_thresh) {
				ggsea_save_overlap_kmer(ctx, overlap_score, t.gid_pos, m.gid_pos_arr[i]);
				continue;
			}

			/* extension */
			ggsea_extend(ctx, t.gid_pos, m.gid_pos_arr[i]);
		}
	}

	/* cleanup iterator */
	gref_iter_clean(iter);

	#if 0
	/* recall saved seeds */
	while(ggsea_finish(ctx) == 0) {
		ggsea_extend(ctx, t.gid_pos, m);
	}
	#endif

	/* traceback results */
	return(ggsea_refine_results(ctx));
}

/**
 * @fn ggsea_aln_free
 */
void ggsea_aln_free(
	ggsea_result_t aln)
{
	return;
}


/* unittests */
/* global configuration */
unittest_config(
	.name = "ggsea",
	.depends_on = { "hmap", "gref", "gaba", "psort" }
);

#define _str(x)		x, strlen(x)
#define _seq(x)		(uint8_t const *)(x), strlen(x)

/* create conf */
unittest()
{
	ggsea_conf_t *conf = ggsea_conf_init(GGSEA_PARAMS(
		.seq_a_format = GABA_ASCII,
		.seq_b_format = GABA_ASCII,
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1),
		.xdrop = 10));

	assert(conf != NULL, "%p", conf);

	ggsea_conf_clean(conf);
}

/* create context */
unittest()
{
	ggsea_conf_t *conf = ggsea_conf_init(GGSEA_PARAMS(
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1),
		.xdrop = 10));

	gref_pool_t *pool = gref_init_pool(GREF_PARAMS( .k = 3 ));
	gref_append_segment(pool, _str("seq1"), _seq("ACGTACGTACGTAACCACGTACGTACGT"));
	gref_acv_t *acv = gref_freeze_pool(pool);
	gref_idx_t *idx = gref_build_index(acv);
	assert(idx != NULL, "%p", idx);

	/* ctx_init fails without reference index */
	ggsea_ctx_t *sea = ggsea_ctx_init(conf, NULL);
	assert(sea == NULL, "%p", sea);

	/* with valid reference index */
	sea = ggsea_ctx_init(conf, idx);
	assert(sea != NULL, "%p", sea);

	/* cleanup */
	ggsea_ctx_clean(sea);
	ggsea_conf_clean(conf);
}

/* single linear sequence */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS( .k = 3 ));
	gref_append_segment(pool, _str("seq1"), _seq("ACGTACGTACGTAACCACGTACGTACGT"));


	assert(1);
}


/**
 * end of ggsea.c
 */
