#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;

const int N = 10;

vector<vector<double>> Matriks(N, vector<double>(N + 1));

double Random_Matriks() {
    return rand() % 11;
}

void Tampil(int pilih) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            if (j == N) {
                if (pilih == 2) {
                    cout << Matriks[i][j];
                }
                if (pilih == 3) {
                    cout << " | " << Matriks[i][j] << "\t";
                }
            } else {
                if (pilih == 1 || pilih == 3) {
                    cout << Matriks[i][j] << "\t";
                }
            }
        }
        cout << endl;
    }
}

double Pembulatan(double value) {
    return round(value * 100.0) / 100.0;
}

double Pengali(double a, double b) {
    return a / b;
}

void OBE(int a, int b, int c, int d) {
    double kali = Pengali(Matriks[a][b], Matriks[c][d]);

    for (int i = 0; i < N + 1; i++) {
        double bantu = Matriks[c][i] * kali;
        Matriks[a][i] = Pembulatan(Matriks[a][i] - bantu);
    }
}

void Eliminasi_Gauss() {
    for (int i = 1; i < N; i++) {
        for (int j = 0; j < i; j++) {
            OBE(i, j, j, j);
        }
    }
    cout << endl;
}

void Hasil() {
    vector<double> x(N);

    x[N - 1] = Matriks[N - 1][N] / Matriks[N - 1][N - 1];

    for (int i = N - 2; i >= 0; i--) {
        double sum = Matriks[i][N];

        for (int j = i + 1; j < N; j++) {
            sum -= Matriks[i][j] * x[j];
        }

        x[i] = sum / Matriks[i][i];
    }

    cout << "Hasil dari persamaan diatas adalah : " << endl;
    cout << fixed << setprecision(2);
    for (int i = 0; i < N; i++) {
        cout << "x" << (i + 1) << " = " << abs(Pembulatan(x[i])) << endl;
    }
}

int main() {
    srand(time(0));

    vector<double> b = {1, -1, 1, -1, 1, -1, 1, -1, 1, -5};

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            if (j == N) {
                Matriks[i][N] = b[i];
            } else {
                Matriks[i][j] = Random_Matriks();
            }
        }
    }

    cout << "\nMatriks Augmented Awal \n";
    Tampil(3);

    Eliminasi_Gauss();
    
    cout << "\nMatriks Segitiga Atas \n";
    Tampil(3);

    Hasil();

    return 0;
}