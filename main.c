#include <stdio.h>
#include </usr/local/include/omp.h>
#include <stdlib.h>
#include <math.h>
#include </opt/homebrew/opt/openblas/include/lapacke.h>
#include </opt/homebrew/opt/openblas/include/cblas.h>

#define N 10
#define M 5

// Fonction de projection de la matrice A sur l'espace de dimension M
void project_matrix(double A[N][N], double P[M][M]) {
    int i, j, k;
    double Q[N][M];

    int m = N, lda = N, ldz = N, info;
    double w[N];
    char jobz = 'V', uplo = 'U';

    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, m, A, lda, w);
    if (info > 0) {
        printf("Le calcul des vecteurs propres a échoué.\n");
        exit(1);
    }
    // Stockage des M vecteurs propres les plus importants dans Q
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            Q[j][i] = A[j][i];
        }
    }

  // Projection de A sur l'espace de dimension M
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            P[i][j] = 0;
            for (k = 0; k < N; k++) {
                P[i][j] += Q[k][i] * A[k][j] * Q[k][j];
            }
        }
    }
}

int main() {
  double A[N][N], P[M][M];
  int i, j;

  // Initialisation de la matrice A avec des valeurs aléatoires
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i][j] = (double)rand()/RAND_MAX;
    }
  }

  // Projection de A sur l'espace de dimension M
  project_matrix(A, P);

  // Affichage de la matrice projetée P
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      printf("%lf ", P[i][j]);
    }
    printf("\n");
  }

  return 0;
}
