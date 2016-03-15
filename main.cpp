#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include "math.h"

#define MAX_VALUE 100

using namespace std;

typedef struct matrix {
    int m, n;
    double ** values;
} * mat;

mat matrix_init(int, int);
void free_matrix(mat);

void matrix_print(mat x) {
    for (int i = 0; i < x->m; i++) {
        for (int j = 0; j < x->n; j++) {
            printf("%.3f\t", x->values[i][j]);
        //    cout << x->values[i][j] << "\t";
        }
        printf("\n");
    }
}

mat matrix_mul(mat x, mat y) {
    if (x->n != y->m)
        cerr << "Matrix dimensions incorrect" << endl;
    mat result = matrix_init(x->m, y->n);
    for (int i = 0; i < x->m; i++) {
        for (int j = 0; j < y->n; j++) {
            result->values[i][j] = 0;
            for (int k = 0; k < x->n; k++) {
                result->values[i][j] += x->values[i][k] * y->values[k][j];
            }
        }
    }
    return result;
}

mat matrix_mul(int x, mat y) {
    mat result = matrix_init(y->m, y->n);
    for (int i = 0; i < y->m; y++) {
        for (int j = 0; j < y->n; y++) {
            result->values[i][j] = x*y->values[i][j];
        }
    }
    return result;
}

mat matrix_transpose(mat x) {
    mat result = matrix_init(x->n, x->m);
    for (int i = 0; i < x->m; i++) {
        for (int j = 0; j < x->n; j++) {
            result->values[j][i] = x->values[i][j];
        }
    }
    return result;
}

double vector_norm(mat v) {
    if (v->n != 1)
        printf("calling vector norm on non-vector");

    double result = 0;
    for (int i = 0; i < v->m; i++) {
        result += pow(v->values[i][0], 2);
    }
    return sqrt(result);
}

pair<mat, mat> QR(mat u) {
    mat q = matrix_init(u->m, u->n);
    mat r = matrix_init(u->m, u->n);

    for (int i = 0; i < u->m; i++) {
        // qi = ui
        for (int j = 0; j < u->n; j++) {
            q->values[j][i] = u->values[j][i];
        }

        for (int j = 0; j < i; j++) {
            // r[j,i] = qj^t * ui
            //mat qj = matrix_init(u->m, 1);
            //for (int k = 0; k < u->m; k++) {
            //    qj->values[k][0] = q->values[k][j];
            //}

            //mat ui = matrix_init(u->m, 1);
            //for (int k = 0; k < u->m; k++) {
            //    ui->values[k][0] = q->values[k][i];
            //}

            //mat qjt = matrix_transpose(qj);

            //mat res = matrix_mul(qjt, ui);
            double res = 0;
            for (int k = 0; k < u->m; k++) {
                res += q->values[k][j] * u->values[k][i];
            }
            //r->values[j][i] = res->values[0][0];
            r->values[j][i] = res;
            //qi = qi - r[j, i]qj
            for (int k = 0; k < u->m; k++) {
                q->values[k][i] = q->values[k][i] - r->values[j][i] * q->values[k][j];
            }

            //free_matrix(qj);
            //free_matrix(qjt);
            //free_matrix(ui);
            //free_matrix(res);
        }
        // r[i.i] = ||qi||
        mat qi = matrix_init(q->m, 1);
        for (int j = 0; j < q->m; j++) {
            qi->values[j][0] = q->values[j][i];
        }
        r->values[i][i] = vector_norm(qi);
        free_matrix(qi);
        // qi = qi/r[i,i]
        if (r->values[i][i]) {
            for (int j = 0; j < q->m; j++) {
                q->values[j][i] = q->values[j][i] / r->values[i][i];
            }
        }
    }
    return make_pair(q, r);
}

pair<mat, mat> QR_ite(mat x) {
    mat ak = x;
    mat qk;
    mat eigenvec = matrix_init(x->m, x->n);
    for (int i = 0; i < x->m; i++) {
        eigenvec->values[i][i] = 1;
    }
    for (int i = 1; i < 30; i++) {
        pair<mat, mat> qr = QR(ak);
        qk = qr.first;
        mat temp = eigenvec;
        eigenvec = matrix_mul(temp, qk);
        free_matrix(temp);
        mat rk = qr.second;
        ak = matrix_mul(rk, qk);
        free_matrix(qr.first);
        free_matrix(qr.second);
        free_matrix(rk);
    }
    return make_pair(ak, eigenvec);
}

mat create_from_file(const char* file_name) {
    ifstream input_file;
    input_file.open(file_name);
    string s;
    int m = 0;
    vector<vector<double> > values(5);
    while (getline(input_file, s)) {
        int i = 0;
        while (i < s.length() && s[i] == ' ') i++;
        if (i == s.length()) break;
        int prev = 0;
        for (int i = 0; i <= s.length(); i++) {
            if (i == s.length() || s[i] == ' ') {
                values[m].push_back(stod(s.substr(prev, i)));
                prev = i + 1;
            }
        }
        m++;
    }
    input_file.close();
    mat return_mat = matrix_init(m, values[0].size());
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < values[0].size(); j++) {
            return_mat->values[i][j] = values[i][j];
        }
    }
    return return_mat;
}

mat create_random(int dimensions) {
    srand((unsigned)time(NULL));
    mat x = matrix_init(dimensions, dimensions);
    for (int i = 0; i < x->m; i++) {
        for (int j = 0; j < x->n; j++) {
            x->values[i][j] = ((double) rand() / (double) RAND_MAX) * MAX_VALUE;
        }
    }
    return x;
}

mat matrix_init(int m, int n) {
    mat x = (mat) malloc(sizeof(struct matrix));
    x->m = m;
    x->n = n;
    x->values = (double **) calloc(m, sizeof(double *));
    for (int i = 0; i < m; i++) {
        x->values[i] = (double *) calloc(n, sizeof(double));
    }
    return x;
}

void free_matrix(mat x) {
    for (int i = 0; i < x->m; i++) {
        free(x->values[i]);
    }
    free(x->values);
    free(x);
}

int main(int argc, const char* argv[]) {
    mat input;
    if (argc > 1) {
        input = create_from_file(argv[1]);
    } else {
        input = create_random(5);
    }

    pair<mat, mat> result = QR_ite(input);

    FILE * fp;
    fp = fopen("output.txt", "w+");
    fprintf(fp, "Eigenvalues:\n");
    for (int i = 1; i <= input->m; i++) {
        fprintf(fp, "%d: %f", i, result.first->values[i-1][i-1]);
    }
    fclose(fp);
    return 0;
}