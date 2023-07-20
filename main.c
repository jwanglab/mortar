#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <zlib.h>
#include <math.h>
#include "kseq.h"
#include "kvec.h"
#include "khash.h"
#include <getopt.h>
#include <string.h>
#include "paf.h"
#include "bed.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(str2int, uint32_t);

void version() {
  printf("mortar version 1.0.0\n   _     _\n  |_|  _|_|\n _|_|_|_|_|\n|_|_|_|_|_|\n");
}

void usage() {
  printf("Build a reference-based consensus from a PAF (Minimap2) alignment with --cs flag\n");
  printf("Usage: mortar [options] <ref> <paf>\n");
  printf("Options:\n");
  printf("  Reference FASTA file\n");
  printf("  PAF alignment of reads to this reference, with --cs flag\n");
  printf("  -b: BED file containing multiplex (tiled) PCR primers used\n");
  printf("      If a BED file is included, primers will be identified and trimmed off before consensus\n");
  printf("  -m: minimum depth below which reference sequence is masked with Ns [default: 20]\n");
  printf("  -d: maximum depth to consider (additional will simply be ignored) [default: 1000]\n");
  printf("  -r: maximum read length (longer reads [potentially chimeric] will simply be ignored) [default: 2700]\n");
  printf("  -l: fraction of depth to call deletion (default: 0.7)\n");
  printf("      the maximum fraction in a contiguous deletion must exceed this value to call a deletion\n");
  printf("      all other sites within a deletion must just be the highest allele\n");
  printf("  -c: output coverage per site to the given file\n");
  printf("  -v: version\n");
  printf("  -h: help\n");
}

static char* ALLELES = "ACGTN-";

typedef struct locus {
  char ref_al;
  uint32_t alt_counts[6]; // = {0,0,0,0,0,0}; // A, C, G, T, N, -(del)
  khash_t(str2int) *ins; // = kh_init(str2int);
} locus;

typedef struct ref {
  locus* loci;
  uint32_t l; // length
} ref;

typedef struct tile {
  uint32_t id;
  uint32_t left_st;
  uint32_t left_en;
  uint32_t right_st;
  uint32_t right_en;
} tile;

// deletion fraction at locus i
float del_frac(ref template, int i) {
  return (float)template.loci[i].alt_counts[5] / (template.loci[i].alt_counts[0] + template.loci[i].alt_counts[1] + template.loci[i].alt_counts[2] + template.loci[i].alt_counts[3] + template.loci[i].alt_counts[4] + template.loci[i].alt_counts[5]);
}

int main(int argc, char *argv[]) {
  char* ref_file = NULL;
  char* paf_file = NULL;
  char* bed_file = NULL;
  char* covg_file = NULL;
  int min_depth = 20;
  int max_depth = 1000;
  int read_max = 2700;
  float del_threshold = 0.7;

  int opt, index;
  opterr = 0;
  while ((opt = getopt(argc, argv, "vhm:d:b:r:l:c:")) != -1) {
    switch (opt) {
      case 'm':
        min_depth = atoi(optarg);
        break;
      case 'd':
        max_depth = atoi(optarg);
        break;
      case 'r':
        read_max = atoi(optarg);
        break;
      case 'l':
        del_threshold = atof(optarg);
        break;
      case 'b':
        bed_file = optarg;
        break;
      case 'c':
        covg_file = optarg;
        break;
      case 'v':
        version();
        return 0;
        break;
      case 'h':
        version();
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'm' || optopt == 'd' || optopt == 'b' || optopt == 'r' || optopt == 'l' || optopt == 'c')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        usage();
        return 1;
    }
  }
  if(optind + 2 > argc) {
    fprintf(stderr, "ERROR: reference FASTA and PAF file are required\n");
    usage();
    return 1;
  }
  ref_file = argv[optind];
  paf_file = argv[optind+1];

  /*
  for (int index = optind; index < argc; index++) {
    printf ("Non-option argument %s\n", argv[index]);
  }
  */

  // generic hashing variables
  khint_t bin, bin2; // hash bin (result of kh_put/get)
  int absent;

  // ------ read primer BED file and put primers into kvec "tiles" ------
	kvec_t(tile) tiles;
	kv_init(tiles);
  if(bed_file) {
    bed_file_t bed = bed_init(bed_file);
    bed_line_t* line = NULL;
    uint32_t tile_id = 0, first_st, first_en;
    uint32_t id;
    while((line = bed_read_line(&bed)) != NULL) {
      //fprintf(stderr, "%s\t%d\t%d\t%s\n", line->chrom, line->st, line->en, line->name);
      size_t lr_st = 0;
      size_t id_st = 0;
      for(int i = strlen(line->name)-1; i >= 0; i--) {
        if(lr_st == 0 && line->name[i] == '_') lr_st = i+1;
        else if(line->name[i] == '_') id_st = i+1;
      }
      id = atoi(line->name+id_st);
      //fprintf(stderr, "%s: %d, %s\n", line->name, id, line->name+lr_st);
      if(strcmp(line->name+lr_st, "LEFT") == 0) {
        if(tile_id != 0) {
          fprintf(stderr, "ERROR: improperly paired primers in tile BED at '%s'\n", line->name);
          return 1;
        }
        tile_id = id;
        first_st = line->st;
        first_en = line->en;
      } else if(strcmp(line->name+lr_st, "RIGHT") == 0) {
        if(tile_id == 0) {
          fprintf(stderr, "ERROR: improperly paired primers in tile BED at '%s'\n", line->name);
          return 1;
        }
        tile t = {tile_id, first_st, first_en, line->st, line->en};
        kv_push(tile, tiles, t);
        tile_id = 0;
      }
      free(line);
    }
    bed_close(&bed);
  }

  // ------ read reference FASTA and create ref structs ------
  gzFile f;
  kseq_t *ks;

  f = gzopen(ref_file, "r");
  if(!f) {
    fprintf(stderr, "Cannot open '%s'\n", ref_file);
    return 1;
  }
  khash_t(str2int) *refmap = kh_init(str2int);
  kvec_t(ref) refs;
  kv_init(refs);

  ks = kseq_init(f); 
  while (kseq_read(ks) >= 0) {
    ref r;
    r.l = ks->seq.l;
    r.loci = calloc(r.l, sizeof(locus));
    for(int i = 0; i < r.l; i++) {
      r.loci[i].ref_al = ks->seq.s[i];
      r.loci[i].ins = kh_init(str2int);
    }
    bin = kh_get(str2int, refmap, ks->name.s);
    if(bin != kh_end(refmap)) { // ref with this name already exists - if allowed to continue, this will produce undefined results and likely crash
      fprintf(stderr, "ERROR: reference sequence '%s' exists more than once, please fix it\n", ks->name.s);
      return 1;
    }
    char* ref_name = malloc(sizeof(char)*(ks->name.l+1));
    ref_name[ks->name.l] = '\0';
    strcpy(ref_name, ks->name.s);
	  bin = kh_put(str2int, refmap, ref_name, &absent);
    kh_val(refmap, bin) = kv_size(refs);
    kv_push(ref, refs, r);
  }


  // ------ main loop: process PAF alignments ------

  paf_file_t *p = paf_open(paf_file);
  if(!p) {
    fprintf(stderr, "ERROR: Cannot open '%s'\n", paf_file);
    return 1;
  }
  paf_rec_t r;
  int ret = paf_read(p, &r);
  uint64_t n = 0;
  char* last = malloc(sizeof(char) * 100);
  last[0] = '\0';
  uint16_t last_m = 0;
  uint32_t lines =  0;
  for(; ret == 0; ret = paf_read(p, &r)) {
    lines++;
    if(lines%10000==0)
      fprintf(stderr, ".", lines);
    if(strcmp(last, r.qn) != 0) {
      last_m = 0;
    }
    if(bed_file && r.ql > read_max) { // too-large, potentially chimeric read - ONLY apply if using BED file
      continue;
    }
    // --cs strings are already normalized to the fw strand, so we don't have to worry about the orientations
    if(r.ml >= last_m) { // primary alignment
      strcpy(last, r.qn);
      last_m = r.ml;
      char* c = r.cs;
      //fprintf(stderr, "%s: %s\n", r.qn, c);
      long int val;
      char curr_op = '\0';
      uint32_t t = r.ts;
      //fprintf(stderr, "%s: t %d : %d\n", r.tn, r.ts, r.te);
      // get reference ID for this target name
      bin = kh_get(str2int, refmap, r.tn);
	    absent = (bin == kh_end(refmap));
      if(absent) {
        fprintf(stderr, "ERROR: PAF target sequence '%s' not found in reference FASTA.", r.tn);
        return 1;
      }
      uint32_t ref_idx = kh_value(refmap, bin);
      ref template = kv_A(refs, ref_idx);
      
      // trim primers based on BED, if provided
      // - this will also implicitly resolve chimeric reads by just discarding all portions outside the first overlapped tile
      uint32_t al_st = r.ts;
      uint32_t al_en = r.te; // exclusive
      if(bed_file) {
        int best_tile = 0; // tile with the greatest overlap
        int best_ovl = 0;
        for(int i = 0; i < kv_size(tiles); i++) {
          // NOTE: all primer end coordinates in the BED are *inclusive*
          int ovl = (kv_A(tiles, i).right_st < r.te ? kv_A(tiles, i).right_st : r.te) - (kv_A(tiles, i).left_en > r.ts ? kv_A(tiles, i).left_en : r.ts);
          //fprintf(stderr, "%d - %d ovl tile %d (%d - %d): %d\n", r.ts, r.te, i, kv_A(tiles, i).left_en, kv_A(tiles, i).right_st, ovl);
          //if(ovl > (kv_A(tiles, i).right_st - kv_A(tiles, i).left_en)*0.8) {
          //  ...
          //  break;
          if(ovl > best_ovl) {
            best_ovl = ovl;
            best_tile = i;
          } else if(ovl == 0 && best_ovl > 0) { // if this tile does not overlap and we already saw one that did, we can stop looking
            break;
          }
        }
        // indicate these cutoff points so that only nt aligning between left_en and right_st are considered
        al_st = kv_A(tiles, best_tile).left_en + 1;
        al_en = kv_A(tiles, best_tile).right_st; // exclusive
      }
      //fprintf(stderr, "%s: al %d - %d, tile %d - %d\n", r.qn, r.ts, r.te, al_st, al_en);

      char* ins_st;
      while(*c) { // check if c is not pointing to the end of the string (null)
        if(curr_op == '+' && (c[0] == '*' || c[0] == ':' || c[0] == '+' || c[0] == '-')) {
          if(t >= al_st && t < al_en) {
            char* ins = malloc(sizeof(char) * (c-ins_st+1));
            ins[c-ins_st] = '\0';
            strncpy(ins, ins_st, c-ins_st);
            //fprintf(stderr, "ref %d (len %d), pos %d (%c), INS: %s\n", ref_idx, template.l, t, template.loci[t].ref_al, ins);
            bin = kh_put(str2int, template.loci[t].ins, ins, &absent); // the problem is NOT the ins string itself - errors in the same place no matter what's inserted here
            if(absent) {
              //kh_key(template.loci[t].ins, bin) = ins; // I don't think this is necessary, the pointer has already been put in the hash if it was absent
              kh_val(template.loci[t].ins, bin) = 1;
            }
            else {
              kh_val(template.loci[t].ins, bin)++;
              free(ins);
            }
          }
          curr_op = '\0';
        }
        if(c[0] == '*') { // mismatch
          c = c + 2;
          if(t >= al_st && t < al_en) {
            int a = 4;
            if(c[0] == 'a') a = 0;
            else if(c[0] == 'c') a = 1;
            else if(c[0] == 'g') a = 2;
            else if(c[0] == 't') a = 3;
            template.loci[t].alt_counts[a]++; // increment alt allele count
          }
          curr_op = '*';
          t++;
        } else if(c[0] == ':') { // match
          val = strtol(c+1, &c, 10);
          //fprintf(stderr, "MATCH %d: %d-%d (loci: %d)\n", val, t, t+val-1, template.l);
          if(t >= al_st && t < al_en) {
            for(int i = t; i < t+val; i++) {
              int a = 4;
              if(template.loci[i].ref_al == 'A') a = 0;
              else if(template.loci[i].ref_al == 'C') a = 1;
              else if(template.loci[i].ref_al == 'G') a = 2;
              else if(template.loci[i].ref_al == 'T') a = 3;
              template.loci[i].alt_counts[a]++; // increment ref allele count
            }
          }
          curr_op = ':';
          t = t + val;
          // we already advanced c here past the match, so we don't need to finish the loop
          continue;
        } else if(c[0] == '-') {
          curr_op = '-';
        } else if(c[0] == '+') {
          curr_op = '+';
          ins_st = c+1;
        } else if(curr_op == '-') {
          //fprintf(stderr, "DEL: %d (loci: %d)\n", t, template.l);
          if(t >= al_st && t < al_en) {
            template.loci[t].alt_counts[5]++; // increment -(del) count
          }
          t++;
        }
        /*
        if(c[0] != '+' && c[0] != '-' && c[0] != ':' && c[0] != '*' && c[0] != 'a' && c[0] != 'g' && c[0] != 'c' && c[0] != 't' && (c[0] < 48 || c[0] > 57)) {
          fprintf(stderr, "UNKNOWN char: %c (%d)\n", c[0], c[0]);
        }
        */

        // if curr_op is '+' we don't need to do anything right now
        // there shouldn't be any more possible cases...
        c++;
      }
    }
  }
  fprintf(stderr, "\n");
  free(last);
  paf_close(p);

  FILE *covg_f = NULL;
  if(covg_file) {
    covg_f = fopen(covg_file,"w");
  }

  // ------ build consensus from allele/ins counts ------
  for (bin = kh_begin(refmap); bin != kh_end(refmap); ++bin) {
		if (kh_exist(refmap, bin)) {
      uint32_t ref_idx = kh_value(refmap, bin);
      ref template = kv_A(refs, ref_idx);
      uint32_t masked = 0;
      fprintf(stdout, ">%s\n", kh_key(refmap, bin));
      fprintf(stderr, "Building consensus of %s\n", kh_key(refmap, bin));

      int gap_start = -1;
      // for each ref seq...:
      for(int i = 0; i < template.l; i++) {
        uint32_t nondel_depth = template.loci[i].alt_counts[0] + template.loci[i].alt_counts[1] + template.loci[i].alt_counts[2] + template.loci[i].alt_counts[3] + template.loci[i].alt_counts[4];
        uint32_t depth = nondel_depth + template.loci[i].alt_counts[5];
        //fprintf(stderr, "%d: %u\n", i, depth);

        if(covg_file) {
          fprintf(covg_f, "%d\n", depth);
        }

        if(depth < min_depth) {
          fprintf(stdout, "N");
          if(gap_start == -1) {
            gap_start = i;
          }
          masked++;
          continue;
        } else {
          if(gap_start != -1) {
            fprintf(stderr, "Gap: [%d, %d)\n", gap_start, i);
            gap_start = -1;
          }
        }

        // insertions:
        uint32_t tot_ins = 0;
        for (bin2 = kh_begin(template.loci[i].ins); bin2 != kh_end(template.loci[i].ins); ++bin2) {
          if (kh_exist(template.loci[i].ins, bin2)) {
            const char* ins = kh_key(template.loci[i].ins, bin2);
            uint32_t ct = kh_val(template.loci[i].ins, bin2);
            if(ct > depth/2) {
              fprintf(stderr, "%d (covg %d): ", i, depth);
              for(int j = (i < 10 ? 0 : i-10); j < i; j++) fprintf(stderr, "%c", template.loci[j].ref_al);
              fprintf(stderr, "[+%s]", ins);
              for(int j = i; j < (i+10 < template.l ? i+10 : template.l); j++) fprintf(stderr, "%c", template.loci[j].ref_al);
              fprintf(stderr, "\n");
              fprintf(stdout, "%s", ins);
            }
            tot_ins = tot_ins + ct;
          }
        }

        //fprintf(stderr, "pos %d: als: %d, %d, %d, %d, %d, %d\n", i, template.loci[i].alt_counts[0], template.loci[i].alt_counts[1], template.loci[i].alt_counts[2], template.loci[i].alt_counts[3], template.loci[i].alt_counts[4], template.loci[i].alt_counts[5]);
        int alt_max = 4;
        int non_del_max = 4;
        for(int j = 0; j < 6; j++) {
          if(template.loci[i].alt_counts[j] > template.loci[i].alt_counts[alt_max]) {
            alt_max = j;
            if(j != 5) { // highest allele excluding deletion
              non_del_max = j;
            }
          }
        }
        if(alt_max == 5) { // deletion
          // count non-deleted homopolymers adjacent and including the target site (l)
          int hp = 0; // homopolymer count
          // consider this allele if it matches the one before or after, otherwise take the one after
          char al = 'N';
          if(i == 0 || i == template.l-1)
            al = template.loci[i].ref_al;
          else
            al = template.loci[i].ref_al == template.loci[i-1].ref_al || template.loci[i].ref_al == template.loci[i+1].ref_al ? template.loci[i].ref_al : template.loci[i+1].ref_al;
          float max_df = del_frac(template, i); // maximum deletion fraction within this deletion-homopolymer block
          if(template.loci[i].ref_al == al && max_df > 0.5) hp++;
          for(int s = i-1; s >= 0 && (del_frac(template, s) > 0.5 || template.loci[s].ref_al == al); s--) {
            float df = del_frac(template, s);
            if(df < 0.5) hp++;
            if(df > max_df) max_df = df;
          }
          for(int s = i+1; s < template.l && (del_frac(template, s) > 0.5 || template.loci[s].ref_al == al); s++) {
            float df = del_frac(template, s);
            if(df < 0.5) hp++;
            if(df > max_df) max_df = df;
          }

          if(hp >= 5 || max_df < del_threshold) {
            fprintf(stdout, "%c", ALLELES[non_del_max]);
            //fprintf(stderr, " [IGNORED]");
          } else {
            // debugging output
            fprintf(stderr, "%d: als (%dX): %d, %d, %d, %d, %d, %d (del %: %.2f, max_df: %.2f), homopolymers: %d -- ", i, depth, template.loci[i].alt_counts[0], template.loci[i].alt_counts[1], template.loci[i].alt_counts[2], template.loci[i].alt_counts[3], template.loci[i].alt_counts[4], template.loci[i].alt_counts[5], (float)template.loci[i].alt_counts[5]/depth*100, max_df, hp);
            for(int j = (i < 10 ? 0 : i-10); j < i; j++) fprintf(stderr, "%c", template.loci[j].ref_al);
            fprintf(stderr, "[-%c]", template.loci[i].ref_al);
            for(int j = i+1; j < (i+10 < template.l ? i+10 : template.l); j++) fprintf(stderr, "%c", template.loci[j].ref_al);
            fprintf(stderr, "\n");
          }
        } else {
          fprintf(stdout, "%c", ALLELES[alt_max]);
          if(ALLELES[alt_max] != template.loci[i].ref_al) {
            fprintf(stderr, "%d\t%c\t%c", i, template.loci[i].ref_al, ALLELES[alt_max]);
            // mark likely APOBEC3 deamination
            if((i > 0 && template.loci[i].ref_al == 'C' && template.loci[i-1].ref_al == 'T' && ALLELES[alt_max] == 'T') || (i < template.l && template.loci[i].ref_al == 'G' && template.loci[i+1].ref_al == 'A' && ALLELES[alt_max] == 'A')) {
              fprintf(stderr, " [APOBEC3 deamination-like]");
            }
            else if((i > 0 && template.loci[i].ref_al == 'C' && template.loci[i-1].ref_al == 'C' && ALLELES[alt_max] == 'T') || (i < template.l && template.loci[i].ref_al == 'G' && template.loci[i+1].ref_al == 'G' && ALLELES[alt_max] == 'A')) {
              fprintf(stderr, " [APOBEC3G deamination-like]");
            }
            fprintf(stderr, "\n");
          }
        }
      }
      fprintf(stdout, "\n");
      fprintf(stderr, "%d sites masked (Ns) - %.2f%% of %u loci\n", masked, ((float)masked*100/template.l), template.l);
    }
  }

  if(covg_file) {
    fclose(covg_f);
  }

  // --- clean up memory ---
	kv_destroy(tiles);

  for (bin = kh_begin(refmap); bin != kh_end(refmap); ++bin) {
		if (kh_exist(refmap, bin)) {
      free((char*)kh_key(refmap, bin));
    }
  }
  kh_destroy(str2int, refmap);

  for(int i = 0; i < kv_size(refs); i++) {
    for(int j = 0; j < kv_A(refs, i).l; j++) {
      for (bin = kh_begin(kv_A(refs, i).loci[j].ins); bin != kh_end(kv_A(refs, i).loci[j].ins); ++bin) {
        if (kh_exist(kv_A(refs, i).loci[j].ins, bin)) {
          free((char*)kh_key(kv_A(refs, i).loci[j].ins, bin));
        }
      }
      kh_destroy(str2int, kv_A(refs, i).loci[j].ins);
    }
    free(kv_A(refs, i).loci);
  }
  kv_destroy(refs);

  return 0;
}

