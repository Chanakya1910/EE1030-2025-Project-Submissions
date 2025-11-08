#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main() {
    int width, height, channels;
    char image[200];
    printf("enter image name:");
    scanf("%s",&image);
    unsigned char *img = stbi_load(image, &width, &height, &channels, 1);
    if (!img) {
        printf("Error: could not load image.\n");
        return 1;
    }

    printf("Loaded %dx%d grayscale image.\n", width, height);
    int m = height, n = width;

    double *A = malloc(sizeof(double) * m * n);
    for (int i = 0; i < m * n; i++) A[i] = img[i];

    double normA = frob(A, m, n);
    printf("enter no.of values of k:");
    int z;
    scanf("%d",&z);
    int k_values[z];
    printf("enter k values:");
    for (int i = 0; i <z ; i++)
    {
        scanf("%d",&k_values[i]);
    }
    
    int count = sizeof(k_values) / sizeof(k_values[0]);

    for (int i = 0; i < count; i++)
        svd(A, m, n, k_values[i], normA);

    stbi_image_free(img);
    free(A);
     return 0;
}

//finds Av ,project back u into row space of A
void Rmult(int m, int n, const double *A, const double *x, double *y) {
    for (int i = 0; i < m; i++) {
    double t = 0;
    for (int j = 0; j < n; j++) t += A[i * n + j] * x[j];
        y[i] = t;
    }
}

//finds v=Atu ,project  into columnspace of A
void Cmult(int m, int n, const double *A, const double *y, double *x) {
    for (int j = 0; j < n; j++) {
    double t = 0;
                for (int i = 0; i < m; i++) t += A[i * n + j] * y[i];
        x[j] = t;
    }
}

//removing singular to compute next singular
void deflate(int m, int n, double *A, const double *u, const double *v, double s) {
    for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
            A[i * n + j] -= s * u[i] * v[j];
}
//claculating singular
int singular(
    int m, int n, const double *A, double *u, double *v, double *s,
    int max_iter, double tol
) {
    double *temp1 = malloc(sizeof(double) * m);
    double *temp2 = malloc(sizeof(double) * n);
    if (!temp1 || !temp2) return -1;

   
    for (int i = 0; i < m; i++) u[i] = (i) + 1;

    double prev = 0;
    for (int i = 0; i < max_iter; i++) {
        Cmult(m, n, A, u, temp2);
        Rmult(m, n, A, temp2, temp1);

        double nv = norm(m, temp1);
        if (nv == 0) break;

        for (int i = 0; i < m; i++) u[i] = temp1[i] / nv;

        if (fabs(nv - prev) < tol * (nv + 1e-9)) break;
        prev = nv;
    }

    Cmult(m, n, A, u, temp2);
    *s = norm(n, temp2);
    if (*s == 0) return -2;

    for (int j = 0; j < n; j++) v[j] = temp2[j] / (*s);

    free(temp1);
    free(temp2);
    return 0;
}
// truncated SVD


 void svd(double *A, int m, int n, int k, double normA) {
    double *U = calloc(k * m, sizeof(double));
    double *V = calloc(k * n, sizeof(double));
    double *S = calloc(k, sizeof(double));
    double *B = malloc(sizeof(double) * m * n);
    memcpy(B, A, sizeof(double) * m * n);

    for (int r = 0; r < k; r++) {
        double *u = U + r * m;
        double *v = V + r * n;
        double s = 0;

        int rc =singular(m, n, B, u, v, &s, 800, 1e-6);
        if (rc != 0) { k = r; break; }

        S[r] = s;
        deflate(m, n, B, u, v, s);
    }

   //new matrix
    double *R = calloc(m * n, sizeof(double));
    for (int r = 0; r < k; r++) {
        double s = S[r];
        double *u = U + r * m;
        double *v = V + r * n;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                R[i * n + j] += s * u[i] * v[j];
    }

   //error
    double diff = 0;
    for (int i = 0; i < m * n; i++) {
        double d = A[i] - R[i];
        diff += d * d;
    }
    double error = (sqrt(diff) / normA)*100;
  printf("---Frobenius error: %.2f ",error);

    unsigned char *out = malloc(m * n);
    for (int i = 0; i < m * n; i++) {
        int val = (int)(R[i] + 0.5);
        if (val < 0) val = 0;
        if (val > 255) val = 255;
        out[i] = (unsigned char)val;
    }
//new image
    char file[64];
    sprintf(file, "compressed k_%d.jpg", k);
    stbi_write_jpg(file, n, m, 1, out, 90);
    printf("Saved %s\n", file);

    free(U); free(V); free(S); free(B); free(R); free(out);
}

//vector norm
double norm(int n, const double *v) {
 double s = 0;
 for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}
//frobenious norm
double frob(const double *A, int m, int n) {
    double s = 0;
    for (int i = 0; i < m * n; i++) s += A[i] * A[i];
    return sqrt(s);
}
