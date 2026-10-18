/*
 * splitrel.c -- exact split-reliability coefficients for two-terminal graphs.
 *
 * Reads graphs in graph6 format (one per line, n <= 12) from stdin.
 * For every graph and every unordered terminal pair {s,t}, prints
 *   <graph6> <s> <t> <m> N_{n-2} N_{n-1} ... N_m
 * where N_i is the number of spanning subgraphs with exactly i edges that have
 * exactly two connected components, one containing s and the other containing t
 * (Brown & McMullin 2023, Theorem 1.1).
 *
 * Method: for every vertex subset S the polynomial c_S(x) = sum_j c_j(G[S]) x^j
 * (c_j = number of connected spanning subgraphs of G[S] with j edges) is obtained
 * by the standard subset recursion
 *    c_S(x) = (1+x)^{e(S)} - sum_{T subset S, v0 in T, T != S} c_T(x) (1+x)^{e(S\T)}
 * with v0 the lowest vertex of S, and then
 *    SR_{s,t}(x) = sum_{A : s in A, t notin A} c_A(x) c_{V\A}(x).
 *
 * Usage: ./splitrel [--maxonly]  < graphs.g6
 *   --maxonly : instead of printing every pair, print for each m the componentwise
 *               maximum vector. This mode requires a single graph order.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 12
#define MAXM (MAXN*(MAXN-1)/2)
#define MAXS (1<<MAXN)

typedef long long ll;

static ll binom[MAXM+1][MAXM+1];
static int n, m;
static unsigned adj[MAXN];
static int ecount[MAXS];
static ll *cpoly[MAXS];      /* cpoly[S][j], j = 0..ecount[S] */
static int clen[MAXS];

static int parse_g6(const char *s) {
    size_t length = strlen(s);
    if (length == 0) return -1;
    int nn = (unsigned char)s[0] - 63;
    if (nn < 1 || nn > MAXN) return -1;
    int bits_total = nn * (nn - 1) / 2;
    size_t data_chars = (size_t)(bits_total + 5) / 6;
    if (length != 1 + data_chars) return -1;
    for (size_t j = 1; j < length; j++) {
        unsigned char c = (unsigned char)s[j];
        if (c < 63 || c > 126) return -1;
    }
    int padding = (6 - bits_total % 6) % 6;
    if (padding && (((unsigned char)s[length-1] - 63) & ((1 << padding) - 1))) return -1;
    n = nn;
    for (int i = 0; i < n; i++) adj[i] = 0;
    int bit = 0, pos = 1, val = 0, nbits = 0;
    m = 0;
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < j; i++) {
            if (nbits == 0) { val = s[pos++] - 63; nbits = 6; }
            bit = (val >> 5) & 1; val <<= 1; nbits--;
            if (bit) { adj[i] |= 1u << j; adj[j] |= 1u << i; m++; }
        }
    }
    return 0;
}

static void compute_polys(void) {
    int full = (1 << n);
    ecount[0] = 0;
    for (int S = 1; S < full; S++) {
        int v = __builtin_ctz(S);
        int rest = S & ~(1 << v);
        ecount[S] = ecount[rest] + __builtin_popcount(adj[v] & rest);
    }
    for (int S = 1; S < full; S++) {
        int e = ecount[S];
        clen[S] = e + 1;
        ll *c = cpoly[S];
        for (int j = 0; j <= e; j++) c[j] = binom[e][j];
        int v = __builtin_ctz(S);
        int rest = S & ~(1 << v);
        /* enumerate proper submasks sub of rest; T = sub | v */
        if (rest) for (int sub = (rest - 1) & rest; ; sub = (sub - 1) & rest) {
            int T = sub | (1 << v);
            int U = S & ~T;
            int eT = ecount[T], eU = ecount[U];
            ll *cT = cpoly[T];
            for (int i = 0; i <= eT; i++) {
                ll ci = cT[i];
                if (!ci) continue;
                for (int k = 0; k <= eU; k++) c[i + k] -= ci * binom[eU][k];
            }
            if (sub == 0) break;
        }
    }
}

int main(int argc, char **argv) {
    if (argc > 2 || (argc == 2 && strcmp(argv[1], "--maxonly") != 0 && strcmp(argv[1], "--conn") != 0)) {
        fprintf(stderr, "usage: splitrel [--maxonly|--conn] < graphs.g6\n");
        return 2;
    }
    int maxonly = (argc > 1 && strcmp(argv[1], "--maxonly") == 0);
    int connonly = (argc > 1 && strcmp(argv[1], "--conn") == 0);
    for (int a = 0; a <= MAXM; a++) { binom[a][0] = 1; for (int b = 1; b <= a; b++) binom[a][b] = binom[a-1][b-1] + (b <= a-1 ? binom[a-1][b] : 0); }
    for (int S = 0; S < MAXS; S++) cpoly[S] = NULL;
    static char line[4096];
    /* for maxonly mode */
    ll maxvec[MAXM+2][MAXM+1]; int have[MAXM+2];
    memset(have, 0, sizeof have);
    long ngraphs = 0;
    int first_order = 0;
    for (;;) {
        size_t L = 0;
        int ch;
        while ((ch = fgetc(stdin)) != EOF && ch != '\n') {
            if (ch == 0 || L + 1 >= sizeof line) {
                fprintf(stderr, "NUL byte or overlong input line\n"); return 2;
            }
            line[L++] = (char)ch;
        }
        if (ch == EOF && L == 0) break;
        line[L] = 0;
        while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (!L) continue;
        if (parse_g6(line) < 0) { fprintf(stderr, "invalid or unsupported graph6 input\n"); return 2; }
        if (maxonly && first_order && first_order != n) { fprintf(stderr, "mixed graph orders in --maxonly\n"); return 2; }
        first_order = n;
        int full = 1 << n;
        for (int S = 0; S < full; S++) if (!cpoly[S]) {
            cpoly[S] = (ll*)malloc(sizeof(ll) * (MAXM + 1));
            if (!cpoly[S]) { fprintf(stderr, "allocation failed\n"); return 2; }
        }
        compute_polys();
        ngraphs++;
        if (connonly) { printf("%s %d", line, m); for (int j = 0; j <= m; j++) printf(" %lld", cpoly[full - 1][j]); printf("\n"); continue; }
        ll sr[MAXM + 1];
        for (int s = 0; s < n; s++) for (int t = s + 1; t < n; t++) {
            for (int j = 0; j <= m; j++) sr[j] = 0;
            for (int A = 0; A < full; A++) {
                if (!((A >> s) & 1) || ((A >> t) & 1)) continue;
                int B = (full - 1) & ~A;
                ll *ca = cpoly[A], *cb = cpoly[B];
                int ea = ecount[A], eb = ecount[B];
                for (int i = 0; i <= ea; i++) { ll ci = ca[i]; if (!ci) continue; for (int k = 0; k <= eb; k++) sr[i + k] += ci * cb[k]; }
            }
            if (!maxonly) {
                printf("%s %d %d %d", line, s, t, m);
                for (int j = n - 2; j <= m; j++) printf(" %lld", sr[j]);
                printf("\n");
            } else {
                if (!have[m]) { have[m] = 1; for (int j = 0; j <= m; j++) maxvec[m][j] = sr[j]; }
                else for (int j = 0; j <= m; j++) if (sr[j] > maxvec[m][j]) maxvec[m][j] = sr[j];
            }
        }
    }
    if (ferror(stdin)) { fprintf(stderr, "input read error\n"); return 2; }
    if (!ngraphs) { fprintf(stderr, "no graphs supplied\n"); return 2; }
    if (maxonly) {
        for (int mm = 0; mm <= MAXM; mm++) if (have[mm]) {
            printf("MAX %d", mm);
            for (int j = n - 2; j <= mm; j++) printf(" %lld", maxvec[mm][j]);
            printf("\n");
        }
    }
    fprintf(stderr, "processed %ld graphs\n", ngraphs);
    return 0;
}
