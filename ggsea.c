
/**
 * @file ggsea.c
 *
 * @brief Graph-to-Graph Seed-and-Extend Alignment
 *
 * @author Hajime Suzuki
 * @date 2016/4/12
 */

#define UNITTEST_UNIQUE_ID			10
#include "unittest.h"

#include <stdint.h>
#include "ggsea.h"
#include "hmap/hmap.h"
#include "psort/psort.h"
#include "tree/tree.h"
#include "gref/gref.h"
#include "gaba/gaba.h"
#include "kvec.h"
#include "sassert.h"
#include "lmm.h"
#include "log.h"

/* inline directive */
#define _force_inline				inline

/* constants */
#define MARGIN_SEQ_SIZE				( 64 )
#define MARGIN_SEQ_OFFSET			( 16 )
#define MARGIN_SEQ_LEN				( 32 )

/* max and min */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )


/* structs */
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
	int64_t overlap_width;			/* overlap filter width */
	struct ggsea_params_s params;
};

/**
 * @struct ggsea_rep_kmer_cont_s
 */
struct ggsea_rep_kmer_cont_s {
	hmap_header_t header;
	uint64_t vec_size;
	kvec_t(struct gref_gid_pos_s) rv;
	kvec_t(struct gref_gid_pos_s) qv;
};

/**
 * @struct ggsea_region_s
 */
struct ggsea_region_s {
	tree_node_t h;					/* (40) header */
	int64_t len;					/* q section length */
	int64_t sp, ep;					/* start and end p coordinate */
	int64_t depth;
	int64_t score;
};
// _static_assert(sizeof(struct ggsea_region_s) == 96);

/**
 * @struct ggsea_overlap_s
 */
struct ggsea_overlap_s {
	int64_t score;
	struct ggsea_region_s *region;
};

/**
 * @struct ggsea_segq_s
 * @brief segment info container (to push into heapqueue)
 */
struct ggsea_segq_s {
	int64_t psum;
	gaba_fill_t const *fill;
	uint32_t rgid;
	uint32_t qgid;
};
_static_assert(sizeof(struct ggsea_segq_s) == 24);

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

	/* overlap filter */
	tree_t *tree;

	/* dp context */
	gaba_dp_t *dp;
	struct ggsea_node_s *node;			/* node info array */
	kvec_t(struct ggsea_segq_s) segq;	/* segment queue */
	uint8_t *margin;
	struct gref_section_s fw_margin, rv_margin;
};


/* conf init / clean */
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

	restore(p.k, 0);
	restore(p.kmer_cnt_thresh, 100);
	restore(p.overlap_thresh, 3);
	restore(p.popcnt_thresh, 25);
	restore(p.score_thresh, 0);			/* score positive */

	#undef restore

	/* malloc mem */
	struct ggsea_conf_s *conf = (struct ggsea_conf_s *)malloc(
		sizeof(struct ggsea_conf_s));
	if(conf == NULL) {
		return(NULL);
	}
	memset(conf, 0, sizeof(struct ggsea_conf_s));

	/* init alignment context */
	conf->gaba = gaba_init(GABA_PARAMS(
		/*.popcnt = p.popcnt_thresh,*/
		.xdrop = p.xdrop,
		.score_matrix = p.score_matrix));
	if(conf->gaba == NULL) {
		free(conf);
		return(NULL);
	}

	/* store constants */
	conf->rep_kmer_hmap_size = 1024;
	conf->max_rep_vec_size = 128;
	conf->overlap_width = 32;
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


/* context init / clean */
/**
 * @fn ggsea_ctx_clean
 */
void ggsea_ctx_clean(
	ggsea_ctx_t *ctx)
{
	if(ctx != NULL) {
		hmap_clean(ctx->rep_kmer); ctx->rep_kmer = NULL;
		gaba_dp_clean(ctx->dp); ctx->dp = NULL;
		kv_hq_destroy(ctx->segq);
		free(ctx->margin); ctx->margin = NULL;
		free(ctx);
	}
	return;
}

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
		sizeof(struct ggsea_rep_kmer_cont_s),
		HMAP_PARAMS(.hmap_size = conf->rep_kmer_hmap_size));
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
	ctx->dp = gaba_dp_init(conf->gaba, NULL, NULL);
	if(ctx->dp == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	#if 0
	/* init node array */
	uint32_t sec_cnt = gref_get_section_count(ref);
	ctx->node = (struct ggsea_node_s *)malloc(sizeof(struct ggsea_node_s) * sec_cnt);
	if(ctx->node == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	#endif

	/* init overlap tree */
	ctx->tree = tree_init(sizeof(struct ggsea_region_s) - sizeof(tree_node_t), NULL);

	/* init queue */
	kv_hq_init(ctx->segq);
	if(kv_ptr(ctx->segq) == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	debug("init, hq_size(%llu)", kv_hq_size(ctx->segq));

	/* init margin seq */
	ctx->margin = (uint8_t *)malloc(2 * sizeof(uint8_t) * MARGIN_SEQ_SIZE);
	if(ctx->margin == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	/* init forward margin section */
	memset(ctx->margin, 0, MARGIN_SEQ_SIZE);
	ctx->fw_margin = (struct gref_section_s){
		.gid = 0xfffc,
		.len = MARGIN_SEQ_LEN,
		.base = ctx->margin
	};

	/* init reverse margin section */
	memset(&ctx->margin[MARGIN_SEQ_SIZE], 0, MARGIN_SEQ_SIZE);
	ctx->rv_margin = (struct gref_section_s){
		.gid = 0xfffd,
		.len = MARGIN_SEQ_LEN,
		.base = &ctx->margin[MARGIN_SEQ_SIZE]
	};
	return(ctx);

_ggsea_ctx_init_error_handler:;
	if(ctx != NULL) {
		ggsea_ctx_clean(ctx);
	}
	return(NULL);
}

/**
 * @fn ggsea_ctx_flush
 */
static _force_inline
void ggsea_ctx_flush(
	struct ggsea_ctx_s *ctx,
	gref_acv_t const *query)
{
	debug("flush called");

	/* set sequence info */
	ctx->q = query;

	/* flush queues */
	kv_hq_clear(ctx->segq);

	/* flush dp context for the new read */
	debug("rlim(%p), qlim(%p)", gref_get_lim(ctx->r), gref_get_lim(ctx->q));
	gaba_dp_flush(ctx->dp, gref_get_lim(ctx->r), gref_get_lim(ctx->q));

	debug("flushed");
	return;
}


/* repetitive kmer filters */
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
	psort_full((void *)ptr, size, sizeof(struct gref_gid_pos_s), 0);

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
	if(kv_size(c->rv) > c->vec_size) {
		kv_size(c->rv) = ggsea_dedup_rep_kmer(ctx, kv_ptr(c->rv), kv_size(c->rv));
	}

	kv_push(c->qv, qpos);
	if(kv_size(c->qv) > c->vec_size) {
		kv_size(c->qv) = ggsea_dedup_rep_kmer(ctx, kv_ptr(c->qv), kv_size(c->qv));
	}

	/* expand vector if spilled */
	if(kv_size(c->rv) > c->vec_size && kv_size(c->qv) > c->vec_size) {
		c->vec_size *= 2;
	}
	return;
}


/* overlap filters (not implemented yet) */
/**
 * @fn ggsea_calc_q
 * @brief calculate q coordinate for use in the overlap filter
 */
static _force_inline
int64_t ggsea_calc_q(
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos,
	uint32_t offset)
{
	union ggsea_q_u {
		int64_t q;
		struct ggsea_q_elems_s {
			uint32_t pos;
			uint32_t gid;
		} e;
	};

	return(((union ggsea_q_u){
		.e = ((struct ggsea_q_elems_s){
			.pos = 0x80000000 ^ (rpos.pos - qpos.pos - offset),
			.gid = rpos.gid ^ ((qpos.gid<<16) | (qpos.gid>>16))
		})
	}).q);
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
	/* calc p (pos) and q (key) */
	int64_t pos = rpos.pos + qpos.pos + ctx->conf.params.k;
	int64_t key = ggsea_calc_q(rpos, qpos, ctx->conf.overlap_width);

	debug("pos(%lld), key(%lld)", pos, key);

	/* retrieve leftmost overlapped section */
	struct ggsea_region_s *n = (struct ggsea_region_s *)
		tree_search_key_right(ctx->tree, key);

	if(n != NULL) {
		debug("n(%p), n->h.key(%lld), n->len(%lld), n->sp(%lld), n->ep(%lld)", n, n->h.key, n->len, n->sp, n->ep);
	}

	/* itarate over regions */
	int64_t depth = INT64_MAX;
	while(n != NULL && n->h.key < (key + ctx->conf.overlap_width)) {
		/* check if the depth exceeds the threshold */
		if(n->sp < pos && pos < n->ep) {
			depth = MIN2(depth, n->depth);
		}

		/* retrieve the next region */
		n = (struct ggsea_region_s *)tree_right(ctx->tree, (tree_node_t *)n);
		if(n != NULL) {
			debug("n(%p), n->h.key(%lld), n->len(%lld), n->sp(%lld), n->ep(%lld)", n, n->h.key, n->len, n->sp, n->ep);
		}
	}
	return((depth == INT64_MAX) ? 0 : depth);
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
	/* discard (not implemented yet) */
	return;
}

/**
 * @fn ggsea_update_overlap_section
 */
static _force_inline
void ggsea_update_overlap_section(
	struct ggsea_ctx_s *ctx,
	struct ggsea_fill_pair_s pair,
	struct gaba_result_s const *r)
{
	for(int64_t i = 0; i < r->slen; i++) {

		/* calc p (pos) and q (key) */
		int64_t sp = r->sec[i].apos + r->sec[i].bpos;
		int64_t ep = sp + r->sec[i].alen + r->sec[i].blen;
		int64_t key = ggsea_calc_q(
			(struct gref_gid_pos_s){
				.gid = r->sec[i].aid,
				.pos = r->sec[i].apos
			},
			(struct gref_gid_pos_s){
				.gid = r->sec[i].bid,
				.pos = r->sec[i].bpos
			},
			ctx->conf.overlap_width);


		/* search overlapping regions */
		struct ggsea_region_s *n = (struct ggsea_region_s *)
			tree_search_key_right(ctx->tree, key);

		if(n != NULL) {
			debug("n(%p), n->h.key(%lld), n->len(%lld), n->sp(%lld), n->ep(%lld)", n, n->h.key, n->len, n->sp, n->ep);
		}

		/* itarate over regions */
		int64_t hit = 0;
		while(n != NULL && n->h.key < (key + ctx->conf.overlap_width)) {
			if(n->sp < (sp + 16) && (ep - 16) < n->ep) {
				hit++;
				n->depth++;
				n->score = MAX2(n->score, r->score);
			}

			/* retrieve the next region */
			n = (struct ggsea_region_s *)tree_right(ctx->tree, (tree_node_t *)n);
			if(n != NULL) {
				debug("n(%p), n->h.key(%lld), n->len(%lld), n->sp(%lld), n->ep(%lld)", n, n->h.key, n->len, n->sp, n->ep);
			}
		}

		if(hit == 0) {
			/* add new region */
			*(n = (struct ggsea_region_s *)r - 1) = (struct ggsea_region_s){
				.h.zero = 0,
				.h.key = key,
				.len = ctx->conf.overlap_width,
				.sp = sp,
				.ep = ep,
				.depth = 1,
				.score = r->score
			};

			debug("n(%p), n->h.key(%lld), n->len(%lld), n->sp(%lld), n->ep(%lld)", n, n->h.key, n->len, n->sp, n->ep);
			tree_insert(ctx->tree, (tree_node_t *)n);
		}
	}
	return;
}

#if 0
/**
 * @fn ggsea_clean_overlap_filter
 */
static
void ggsea_free_tree_elem(
	tree_node_t *node)
{
	return;
}
static _force_inline
void ggsea_clean_overlap_filter(
	struct ggsea_ctx_s *ctx)
{
	tree_iterate(ctx->tree, )
	return;
}
#endif

/* graph traverse functions */
/**
 * @macro _rup, _qup
 * @brief extract update flag
 */
#define _rup(_fill)			( (_fill)->status & GABA_STATUS_UPDATE_A )
#define _qup(_fill)			( (_fill)->status & GABA_STATUS_UPDATE_B )
#define _rqup(_fill)		( (_fill)->status & (GABA_STATUS_UPDATE_A | GABA_STATUS_UPDATE_B) )
#define _term(_fill)		( (_fill)->status & GABA_STATUS_TERM )

/**
 * @fn ggsea_extend_update_queue_r
 */
static _force_inline
gaba_fill_t const *ggsea_extend_update_queue_r(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	gaba_fill_t const *max,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	struct gref_link_s rlink = gref_get_link(ctx->r, rsec->gid);

	/* if len == 0, push margin */
	if(rlink.len == 0) {
		debug("r reached leaf");
		/* update rsec with fw_margin */
		rsec = &ctx->fw_margin;

		/* fill loop */
		while(1) {
			fill = gaba_dp_fill(ctx->dp, fill,
				(struct gaba_section_s *)rsec,
				(struct gaba_section_s *)qsec);
			debug("status(%x), max(%lld), r(%u), q(%u)",
				fill->status, fill->max, rsec->gid, qsec->gid);
			max = (fill->max > max->max) ? fill : max;
			if(_term(fill) | _rup(fill)) { break; }
			qsec = (_qup(fill) == 0) ? qsec : &ctx->rv_margin;
		}
		debug("r leaf finished, max(%lld)", max->max);
		return(max);
	}

	/* push section pairs */
	for(int64_t i = 0; i < rlink.len; i++) {
		kv_hq_push(ctx->segq, ((struct ggsea_segq_s){
			.psum = (int64_t)fill->psum,
			.fill = fill,
			.rgid = rlink.gid_arr[i],
			.qgid = qsec->gid
		}));
	}
	return(max);
}

/**
 * @fn ggsea_extend_update_queue_q
 */
static _force_inline
gaba_fill_t const *ggsea_extend_update_queue_q(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	gaba_fill_t const *max,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	struct gref_link_s qlink = gref_get_link(ctx->q, qsec->gid);

	/* if len == 0, push margin */
	if(qlink.len == 0) {
		debug("q reached leaf");
		/* update qsec with rv_margin */
		qsec = &ctx->rv_margin;

		/* fill loop */
		while(1) {
			fill = gaba_dp_fill(ctx->dp, fill,
				(struct gaba_section_s *)rsec,
				(struct gaba_section_s *)qsec);
			debug("status(%x), max(%lld), r(%u), q(%u)",
				fill->status, fill->max, rsec->gid, qsec->gid);
			max = (fill->max > max->max) ? fill : max;

			if(_term(fill) | _qup(fill)) { break; }
			rsec = (_rup(fill) == 0) ? rsec : &ctx->fw_margin;
		}
		debug("q leaf finished, max(%lld)", max->max);
		return(max);
	}

	/* push section pair if all the incoming edges are scaned */
	for(int64_t i = 0; i < qlink.len; i++) {
		kv_hq_push(ctx->segq, ((struct ggsea_segq_s){
			.psum = (int64_t)fill->psum,
			.fill = fill,
			.rgid = rsec->gid,
			.qgid = qlink.gid_arr[i]
		}));
	}
	return(max);
}

/**
 * @fn ggsea_extend_update_queue_rq
 */
static _force_inline
gaba_fill_t const *ggsea_extend_update_queue_rq(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	gaba_fill_t const *max,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	/* update sections */
	struct gref_link_s rlink = gref_get_link(ctx->r, rsec->gid);
	struct gref_link_s qlink = gref_get_link(ctx->q, qsec->gid);

	if((rlink.len | qlink.len) == 0) {
		debug("r and q reached leaf");
		/* fill once */
		fill = gaba_dp_fill(ctx->dp, fill,
			(struct gaba_section_s *)&ctx->fw_margin,
			(struct gaba_section_s *)&ctx->rv_margin);
		debug("status(%x), max(%lld), r(%u), q(%u)",
			fill->status, fill->max, rsec->gid, qsec->gid);
		max = (fill->max > max->max) ? fill : max;
		debug("rq leaf finished, max(%lld)", max->max);
		return(max);
	}

	/* push section pair if all the incoming edges are scaned */
	for(int64_t i = 0; i < rlink.len; i++) {
		for(int64_t j = 0; j < qlink.len; j++) {
			kv_hq_push(ctx->segq, ((struct ggsea_segq_s){
				.psum = (int64_t)fill->psum,
				.fill = fill,
				.rgid = rlink.gid_arr[i],
				.qgid = qlink.gid_arr[j]
			}));
		}
	}
	return(max);
}

/**
 * @fn ggsea_extend_update_queue
 */
static _force_inline
gaba_fill_t const *ggsea_extend_update_queue(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	gaba_fill_t const *max,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	uint64_t rup = _rup(fill) != 0;
	uint64_t qup = _qup(fill) != 0;

	debug("extend_update, r(%llu), q(%llu)", rup, qup);

	if(rup && qup) {
		return(ggsea_extend_update_queue_rq(ctx, fill, max, rsec, qsec));
	} else if(rup) {
		return(ggsea_extend_update_queue_r(ctx, fill, max, rsec, qsec));
	} else if(qup) {
		return(ggsea_extend_update_queue_q(ctx, fill, max, rsec, qsec));
	} else {
		/* never reaches here */
	}
	return(max);
}

/**
 * @fn ggsea_extend_intl
 */
static _force_inline
gaba_fill_t const *ggsea_extend_intl(
	struct ggsea_ctx_s *ctx,
	struct gref_section_s const *rsec,
	uint32_t rpos,
	struct gref_section_s const *qsec,
	uint32_t qpos)
{
	debug("seed: r(%u, %u), q(%u, %u)", rsec->gid, rpos, qsec->gid, qpos);

	/* flush queue */
	kv_hq_clear(ctx->segq);

	/* fill the first section */
	gaba_fill_t const *max = NULL;
	gaba_fill_t const *fill = max = gaba_dp_fill_root(ctx->dp,
		(struct gaba_section_s *)rsec, rpos,
		(struct gaba_section_s *)qsec, qpos);

	debug("root: status(%x), max(%lld), r(%u, %u), q(%u, %u)",
		fill->status, fill->max, rsec->gid, rpos, qsec->gid, qpos);

	/* check xdrop term */
	if((fill->status & GABA_STATUS_TERM) != 0) {
		return(max);
	}

	/* update first joint */
	max = ggsea_extend_update_queue(ctx, fill, max, rsec, qsec);

	/* loop */
	while(kv_hq_size(ctx->segq) > 0) {
		struct ggsea_segq_s seg = kv_hq_pop(ctx->segq);
		if(seg.fill == NULL) {
			debug("fill == NULL (unexpected NULL pointer detected)");
			break;
		}

		/**
		 * Lazy merge detection comes here.
		 * tails at the head of the queue is scaned and overlapping pairs are
		 * detected in the lazy strategy. the detected pairs will be merged
		 * by ggsea_merge_tails function.
		 */

		/* extend */
		gaba_fill_t *fill = gaba_dp_fill(ctx->dp, seg.fill,
			(struct gaba_section_s const *)gref_get_section(ctx->r, seg.rgid),
			(struct gaba_section_s const *)gref_get_section(ctx->q, seg.qgid));

		/* update max */
		debug("check max, max(%lld), prev_max(%lld)", fill->max, max->max);
		max = (fill->max > max->max) ? fill : max;

		/* check xdrop term */
		if((fill->status & GABA_STATUS_TERM) == 0) {
			max = ggsea_extend_update_queue(ctx, fill, max, rsec, qsec);
			debug("queue updated, max(%lld)", max->max);
		}
	}

	debug("extend finished, max(%lld)", max->max);
	return(max);
}

/**
 * @fn ggsea_extend
 */
static _force_inline
struct ggsea_fill_pair_s ggsea_extend(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* forward section */
	debug("forward extend");
	struct gref_section_s const *rfsec = gref_get_section(ctx->r, rpos.gid);
	struct gref_section_s const *qfsec = gref_get_section(ctx->q, qpos.gid);
	gaba_fill_t const *fw_max = ggsea_extend_intl(ctx,
		rfsec, rpos.pos,
		qfsec, qpos.pos);

	/* reverse section */
	debug("reverse extend");
	struct gref_section_s const *rrsec = gref_get_section(ctx->r, gref_rev_gid(rpos.gid));
	struct gref_section_s const *qrsec = gref_get_section(ctx->q, gref_rev_gid(qpos.gid));
	gaba_fill_t const *rv_max = ggsea_extend_intl(ctx,
		rrsec, rrsec->len - rpos.pos,
		qrsec, qrsec->len - qpos.pos);

	debug("fw_max(%lld), rv_max(%lld), max(%lld)",
		fw_max->max, rv_max->max, fw_max->max + rv_max->max);
	/* return max pair */
	return((struct ggsea_fill_pair_s){
		.fw = fw_max,
		.rv = rv_max
	});
}


/**
 * @fn ggsea_align
 */
struct ggsea_result_s ggsea_align(
	ggsea_ctx_t *_ctx,
	gref_acv_t const *query)
{
	struct ggsea_ctx_s *ctx = (struct ggsea_ctx_s *)_ctx;

	debug("align entry");

	/* flush ctxing buffer */
	ggsea_ctx_flush(ctx, query);

	/* seeding */
	gref_iter_t *iter = gref_iter_init(ctx->q);

	debug("iter inited");


	kvec_t(struct gaba_result_s const *) aln;
	kv_init(aln);

	struct gref_kmer_tuple_s t;
	while((t = gref_iter_next(iter)).kmer != GREF_ITER_KMER_TERM) {
		
		/* match */
		struct gref_match_res_s m = gref_match_2bitpacked(ctx->r, t.kmer);

		debug("kmer(%llx), gid(%u), pos(%u), m.len(%lld)",
			t.kmer, t.gid_pos.gid, t.gid_pos.pos, m.len);

		/* check if kmer is repetitive */
		if(m.len > ctx->conf.params.kmer_cnt_thresh) {
			ggsea_save_rep_kmer(ctx, t.kmer, m.gid_pos_arr[0], t.gid_pos);
			continue;
		}

		/* for all positions on ref matched with the kmer */
		for(int64_t i = 0; i < m.len; i++) {
			/* apply overlap filter */
			int64_t score = ggsea_overlap_filter(ctx, m.gid_pos_arr[i], t.gid_pos);

			debug("overlap(%lld, %lld)", score, ctx->conf.params.overlap_thresh);

			if(score >= ctx->conf.params.overlap_thresh) {
				ggsea_save_overlap_kmer(ctx, score, m.gid_pos_arr[i], t.gid_pos);
				continue;
			}

			/* extension */
			struct ggsea_fill_pair_s pair = ggsea_extend(ctx, m.gid_pos_arr[i], t.gid_pos);
			if(pair.fw->max + pair.rv->max <= ctx->conf.params.score_thresh) {
				continue;
			}

			/* traceback */
			struct gaba_result_s const *r = gaba_dp_trace(ctx->dp, pair.fw, pair.rv, NULL);
			kv_push(aln, r);

			/* update overlap filter */
			ggsea_update_overlap_section(ctx, pair, r);
		}
	}

	/* cleanup iterator */
	gref_iter_clean(iter);
	debug("done. %llu alignments generated", kv_size(aln));

	return((struct ggsea_result_s){
		.ref = ctx->r,
		.query = ctx->q,
		.aln = kv_ptr(aln),
		.cnt = kv_size(aln)
	});
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
	.depends_on = { "hmap", "psort", "tree", "gref", "gaba" }
);

#define _str(x)		x, strlen(x)
#define _seq(x)		(uint8_t const *)(x), strlen(x)

#if 0
/* create conf */
unittest()
{
	ggsea_conf_t *conf = ggsea_conf_init(GGSEA_PARAMS(
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
#endif
/* omajinais */
#define with_default_conf() \
	.init = (void *(*)(void *))ggsea_conf_init, \
	.clean = (void (*)(void *))ggsea_conf_clean, \
	.params = (void *)GGSEA_PARAMS( \
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1), \
		.xdrop = 10)

#define omajinai() \
	ggsea_conf_t *conf = (ggsea_conf_t *)ctx;
#if 0
/* omajinai test */
unittest(with_default_conf())
{
	omajinai();
	assert(conf != NULL, "%p", conf);
}

/* single linear sequence */
unittest(with_default_conf())
{
	omajinai();

	/* build sequence pool */
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS( .k = 3 ));
	gref_append_segment(pool, _str("seq1"), _seq("ACGTACGTACGTAACCACGTACGTACGT"));
	gref_idx_t *idx = gref_build_index(gref_freeze_pool(pool));

	/* build ggsea context */
	ggsea_ctx_t *sea = ggsea_ctx_init(conf, idx);

	/* align */
	ggsea_result_t r = ggsea_align(sea, (gref_acv_t *)idx);
	assert(r.aln != NULL, "%p", r.aln);

	ggsea_aln_free(r);
	ggsea_ctx_clean(sea);
	gref_clean(idx);
}
#endif
/* longer sequence */
unittest(with_default_conf())
{
	omajinai();

	/* build reference index object */
	gref_pool_t *rpool = gref_init_pool(GREF_PARAMS( .k = 14 ));
	gref_append_segment(rpool, _str("ref1"),
		_seq("CTCACCTCGCTCAAAAGGGCTGCCTCCGAGCGTGTGGGCGAGGACAACCGCCCCACAGTCAAGCTCGAATGGGTGCTATTGCGTAGCTAGGACCGGCACT"));
	gref_idx_t *ref = gref_build_index(gref_freeze_pool(rpool));

	/* build query sequence iterator */
	gref_pool_t *qpool = gref_init_pool(GREF_PARAMS( .k = 14 ));
	gref_append_segment(qpool, _str("query1"),
		_seq("GGCTGCCTCCGAGCGTGTGGGCGAGGACAACCGCCCCACAGTCAAGCTCGAA"));
	gref_acv_t *query = gref_freeze_pool(qpool);


	/* build ggsea context */
	ggsea_ctx_t *sea = ggsea_ctx_init(conf, ref);

	/* align */
	ggsea_result_t r = ggsea_align(sea, query);
	assert(r.aln != NULL, "%p", r.aln);

	/* cleanup */
	ggsea_aln_free(r);
	ggsea_ctx_clean(sea);
	gref_clean(ref);
	gref_clean(query);
}

/**
 * end of ggsea.c
 */
