/*
 * https://github.com/lh3/miniasm/blob/master/paf.h
 * commit ff3f28a
 */

#ifndef PAF_PAF_H
#define PAF_PAF_H

#include <stdint.h>
#include <sys/types.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
  size_t l, m;
  char *s;
} kstring_t;
#endif

typedef struct {
  void *fp;
  kstring_t buf;
} paf_file_t;

//cluster_0_11	1504	0	1504	-	NZ_CP058220.1	1542	31	1535	1504	1504	0	NM:i:0	ms:i:3008	AS:i:3008	nn:i:0	tp:A:P	cm:i:178	s1:i:1081	s2:i:1081	de:f:0	rl:i:672	cg:Z:1504M	cs:Z::1504
typedef struct {
  const char *qn, *tn; // these point to the input string; NOT allocated
  uint32_t ql, qs, qe, tl, ts, te;
  uint32_t ml:31, rev:1, bl;
  uint8_t mq;
  const char *cigar; // also points to the input string; NOT allocated
  char *cs; // also points to the input string; NOT allocated
} paf_rec_t;

#ifdef __cplusplus
extern "C" {
#endif

paf_file_t *paf_open(const char *fn);
int paf_close(paf_file_t *pf);
int paf_read(paf_file_t *pf, paf_rec_t *r);

#ifdef __cplusplus
}
#endif

#endif
