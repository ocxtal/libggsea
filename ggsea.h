
/**
 * @file ggsea.h
 *
 * @brief Graph-to-Graph Seed-and-Extend Alignment
 *
 * @author Hajime Suzuki
 * @date 2016/4/12
 */
#ifndef _GGSEA_H_INCLUDED
#define _GGSEA_H_INCLUDED

#include <stdint.h>
#include "gref/gref.h"
#include "gaba/gaba.h"


/* types */
/**
 * @type ggsea_conf_t
 */
typedef struct ggsea_conf_s ggsea_conf_t;

/**
 * @type ggsea_ctx_t
 */
typedef struct ggsea_ctx_s ggsea_ctx_t;

/**
 * @struct ggsea_params_s
 */
struct ggsea_params_s {
	/** score parameters */
	int16_t xdrop;
	gaba_score_t const *score_matrix;

	/* repetitive kmer filter */
	int64_t kmer_cnt_thresh;		/* kmer count threshold */

	/* popcnt filter thresh */
	int64_t popcnt_thresh;			/* first threshold */
	int64_t popcnt_low_thresh;		/* rescue threshold */

	/* overlap filter thresh */
	int64_t overlap_thresh;			/* depth */

	/* number of threads in parallel radix sort */
	uint32_t num_threads;
};
typedef struct ggsea_params_s ggsea_params_t;

/**
 * @macro GGSEA_PARAMS
 * @brief utility macro for gaba_init, see example on header.
 */
#define GGSEA_PARAMS(...)			( &((struct ggsea_params_s const) { __VA_ARGS__ }) )

/**
 * @struct ggsea_result_s
 */
struct ggsea_result_s {
	struct gaba_result_s const *const *aln;
	int64_t cnt;
};
typedef struct ggsea_result_s ggsea_result_t;


/* functions */

/**
 * @fn ggsea_conf_init
 * @brief create configuration object
 */
ggsea_conf_t *ggsea_conf_init(
	ggsea_params_t const *params);

/**
 * @fn ggsea_conf_clean
 * @brief cleanup configuration object
 */
void ggsea_conf_clean(
	ggsea_conf_t *conf);



#endif /* _GGSEA_H_INCLUDED */
/**
 * end of ggsea.h
 */
