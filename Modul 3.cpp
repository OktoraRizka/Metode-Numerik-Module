#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// batas nol kecil
const double DELTA = 1e-12;

// Menampilkan isi array
void printArray(double arr[], int len) {
    for (int i = 0; i < len; i++) {
        cout << setw(15) << fixed << setprecision(7) << arr[i];
    }
    cout << endl;
}

// Hitung dot product untuk keperluan back substitution
double dotPart(double row[], double sol[], int from, int to) {
    double total = 0.0;
    for (int i = from; i < to; i++) {
        total += row[i] * sol[i];
    }
    return total;
}

// Eliminasi Gauss (versi logika baru)
double* solveGauss(double M[][3], double V[], int n) {

    for (int col = 0; col < n; col++) {

        // cari pivot maksimum
        int pivot = col;
        for (int r = col + 1; r < n; r++) {
            if (fabs(M[r][col]) > fabs(M[pivot][col])) {
                pivot = r;
            }
        }

        // tukar baris bila diperlukan
        for (int c = 0; c < n; c++) {
            swap(M[col][c], M[pivot][c]);
        }
        swap(V[col], V[pivot]);

        if (fabs(M[col][col]) < DELTA) {
            throw runtime_error("Matrix tidak dapat diselesaikan (singular).");
        }

        // eliminasi baris bawahnya
        for (int r = col + 1; r < n; r++) {
            double factor = M[r][col] / M[col][col];
            V[r] -= factor * V[col];
            for (int c = col; c < n; c++) {
                M[r][c] -= factor * M[col][c];
            }
        }
    }

    // back substitution
    double* sol = new double[n];
    for (int r = n - 1; r >= 0; r--) {
        sol[r] = (V[r] - dotPart(M[r], sol, r + 1, n)) / M[r][r];
    }
    return sol;
}

// menghitung matriks koefisien regresi
void buildCoeff(double pts[][2], int deg, double out[][3], int N) {
    for (int r = 0; r < deg; r++) {
        for (int c = 0; c < deg; c++) {
            double sumX = 0;
            for (int i = 0; i < N; i++) {
                sumX += pow(pts[i][0], r + c);
            }
            out[r][c] = sumX;
        }
    }
    out[0][0] = N;
}

// menghitung konstanta regresi
void buildConst(double pts[][2], int deg, double C[], int N) {
    for (int i = 0; i < deg; i++) {
        double total = 0;
        for (int j = 0; j < N; j++) {
            total += pow(pts[j][0], i) * pts[j][1];
        }
        C[i] = total;
    }
}

// interpolasi lagrange
double lagrangeEval(double X[], double Y[], double xp, int N) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        double L = 1.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                L *= (xp - X[j]) / (X[i] - X[j]);
            }
        }
        res += L * Y[i];
    }
    return res;
}

// tabel lagrange
void showLagrangeTable(double X[], double Y[], int N, double xq, double val) {
    cout << "+------------+--------------------+--------------------+\n";
    cout << "|    X[i]    |        Y[i]        |   Hasil Lagrange   |\n";
    cout << "+------------+--------------------+--------------------+\n";

    for (int i = 0; i < N; i++) {
        cout << "| " << setw(10) << X[i] 
             << " | " << setw(18) << Y[i] 
             << " | " << setw(18) << val 
             << " |\n";
    }

    cout << "+------------+--------------------+--------------------+\n";
    cout << "Nilai Lagrange di z = " << xq << " adalah " << val << endl;
}

int main() {

    double sampel[][2] = {
        {1,1.5577},{2,1.2131},{3,0.9447},{4,0.7358},
        {5,0.5730},{6,0.4462},{7,0.3476},{8,0.2706}
    };

    int N = sizeof(sampel)/sizeof(sampel[0]);

    double* solusi;
    double err = 0;

    // ===================== REGRESI LINIER ==========================
    int deg = 2;
    double A[2][3], B[2];

    buildCoeff(sampel, deg, A, N);
    buildConst(sampel, deg, B, N);

    cout << "\n=== Regresi Linier (a + bx) ===\n";
    printArray(B, deg);

    solusi = solveGauss(A, B, deg);
    cout << "Persamaan: " << solusi[1] << " + " << solusi[0] << "x\n";

    err = 0;
    for (int i = 0; i < N; i++) {
        err += pow((solusi[1] + solusi[0]*sampel[i][0] - sampel[i][1]), 2);
    }
    cout << "Error: " << sqrt(err/N) << endl;

    delete[] solusi;

    // ===================== REGRESI KUADRAT =========================
    deg = 3;
    double A2[3][3], B2[3];

    buildCoeff(sampel, deg, A2, N);
    buildConst(sampel, deg, B2, N);

    cout << "\n=== Regresi Kuadrat (ax^2 + bx + c) ===\n";
    printArray(B2, deg);
    solusi = solveGauss(A2, B2, deg);

    cout << "Persamaan: " << solusi[2] << "x^2 + " 
         << solusi[1] << "x + " << solusi[0] << endl;

    err = 0;
    for (int i = 0; i < N; i++) {
        err += pow((solusi[2]*pow(sampel[i][0],2) + solusi[1]*sampel[i][0] + solusi[0] - sampel[i][1]), 2);
    }
    cout << "Error: " << sqrt(err/N) << endl;

    delete[] solusi;

    // ===================== REGRESI EKSPONENSIAL ====================
    double sampelLog[8][2];
    for (int i = 0; i < N; i++) {
        sampelLog[i][0] = sampel[i][0];
        sampelLog[i][1] = log(sampel[i][1]);
    }

    deg = 2;
    double A3[2][3], B3[2];

    buildCoeff(sampelLog, deg, A3, N);
    buildConst(sampelLog, deg, B3, N);

    cout << "\n=== Regresi Eksponensial (ae^(bx)) ===\n";
    printArray(B3, deg);

    solusi = solveGauss(A3, B3, deg);

    cout << "Persamaan: " << exp(solusi[0]) << " * e^(" << solusi[1] << "x)\n";
    delete[] solusi;

    // ===================== LAGRANGE ================================
    double Xval[] = {1,2,3,4,5,6,7,8};
    double Yval[] = {1.5577,1.2131,0.9447,0.7358,0.5730,0.4462,0.3476,0.2706};
    double z = 5.5;

    cout << "\n=== Interpolasi Lagrange ===\n";
    double Lout = lagrangeEval(Xval, Yval, z, 8);

    showLagrangeTable(Xval, Yval, 8, z, Lout);

    return 0;
}
