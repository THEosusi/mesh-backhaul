double dot(const std::vector<double>& A, const std::vector<double>& B);

void matVecMult(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& result);

int BiCGSTAB(const std::vector<std::vector<double>>& press_coeff, const std::vector<double>& data_all, std::vector<double>& output, int max_iter, double eps);

void BiCGSTABRes(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, std::vector<double>& r, size_t n);

double dot(const std::vector<double>& A, const std::vector<double>& B) {
    double result = 0.0;
    for (size_t i = 0; i < A.size(); ++i) {
        result += A[i] * B[i];
    }
    return result;
}

void matVecMult(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& result) {
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i) {
        result[i] = dot(A[i], x);
    }
}

int BiCGSTAB(const std::vector<std::vector<double>>& press_coeff, const std::vector<double>& data_all, std::vector<double>& output, int node, int max_iter, double eps) {
    size_t n = node;
    std::vector<double> r(n), p(n), v(n), s(n), t(n);
    output.assign(n, 0.0);

    BiCGSTABRes(press_coeff, output, data_all, r, n);
    std::vector<double> r0 = r;

    double rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    double res0 = sqrt(dot(r, r));
    if (res0 < eps) {
        return 0;
    }

    for (int k = 0; k < max_iter; ++k) {
        double rho1 = dot(r0, r);
        double beta = (rho1 / rho) * (alpha / omega);
        rho = rho1;

        for (size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        matVecMult(press_coeff, p, v);
        alpha = rho / dot(r0, v);

        for (size_t i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        double res1 = sqrt(dot(s, s));
        if (res1 < eps) {
            for (size_t i = 0; i < n; ++i) {
                output[i] += alpha * p[i];
            }
            return k + 1;
        }

        matVecMult(press_coeff, s, t);
        omega = dot(t, s) / dot(t, t);

        for (size_t i = 0; i < n; ++i) {
            output[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        double res2 = sqrt(dot(r, r));
        if (res2 < eps) {
            return k + 1;
        }
    }

    return -1;
}

void BiCGSTABRes(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, std::vector<double>& r, size_t n) {
    std::vector<double> Ax(n);
    matVecMult(A, x, Ax);
    for (size_t i = 0; i < n; ++i) {
        r[i] = b[i] - Ax[i];
    }
}