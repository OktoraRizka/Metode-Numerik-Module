#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace std;

double persamaan(double x) {
    return -pow(x, 3) + pow(x, 2) - 2 * x + (x + 2) * exp(-2 * x);
}

double TitikPotongY(double x) {
    return pow(x, 2) - (2 * x) + 1;
}

void MetodeSecant(double x0, double x1, double error, int max_iterasi) {
    double fx0 = persamaan(x0);
    double fx1 = persamaan(x1);
    double x_baru, fx_baru;
    int iter = 0;

    cout << "=========================================================================================================================================\n";
    cout << setw(10) << " Iterasi " << " | " << setw(15) << "x0" << " | " << setw(15) << "x1" << " |"
         << setw(15) << "h(x0)" << " | " << setw(15) << "h(x1)" << " | " << setw(15) << "x_baru"
         << " | " << setw(15) << "h(x_baru)" << " | " << setw(15) << "error" << " | " << endl;
    cout << "=========================================================================================================================================\n";

    for (iter = 1; iter <= max_iterasi; iter++) {
        if (abs(fx1 - fx0) < 1e-15) {
            cout << "Selisih h(x1) dan h(x0) terlalu kecil, berhenti untuk menghindari pembagian dengan nol." << endl;
            break;
        }

        x_baru = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        fx_baru = persamaan(x_baru);

        cout << fixed << setprecision(7);
        cout << setw(10) << iter << " | " << setw(15) << x0 << " | " << setw(15) << x1 << " | "
             << setw(15) << fx0 << " | " << setw(15) << fx1 << " | " << setw(15) << x_baru
             << " | " << setw(15) << fx_baru << " | " << setw(15) << abs(fx_baru) << " | " << endl;

        x0 = x1;
        fx0 = fx1;
        x1 = x_baru;
        fx1 = fx_baru;

        if (abs(fx1) < error) {
            break;
        }
    }

    double y = TitikPotongY(x1);
    cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;

    if (iter > max_iterasi) {
        cout << "\nMetode tidak konvergen setelah " << max_iterasi << " iterasi." << endl;
    } else {
        // Baris ini sudah ada dan juga informatif
        cout << "\nMetode konvergen setelah " << iter << " iterasi." << endl;
    }

    cout << fixed << setprecision(7);
    cout << "Nilai error terakhir\t\t: " << abs(fx1) << endl;
    cout << "Akar x yang ditemukan\t\t: " << x1 << endl;
    cout << "Nilai y pada titik potong\t: " << y << endl;
    cout << "Berpotongan pada titik (x,y)\t: (" << x1 << ", " << y << ")\n";
    
    // -- BARIS TAMBAHAN UNTUK KETERANGAN SPESIFIK --
    cout << "Titik potong ditemukan pada iterasi ke-" << iter << "." << endl;
}

int main() {
    double x0, x1, error;
    int iterasi;

    cout << "===== METODE SECANT =====\n";
    cout << "Masukkan tebakan awal x0 : "; cin >> x0;
    cout << "Masukkan tebakan awal x1 : "; cin >> x1;
    cout << "Error toleransi          : "; cin >> error;
    cout << "Jumlah iterasi maksimum  : "; cin >> iterasi;

    MetodeSecant(x0, x1, error, iterasi);

    return 0;
}