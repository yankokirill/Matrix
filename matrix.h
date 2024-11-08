#include <algorithm>
#include <array>
#include <complex>
#include <iostream>
#include <vector>

namespace {
    const long double pi = 2 * asinl(1);

    size_t to_pow2(size_t n) {
        size_t m = 1;
        while (m < n) {
            m *= 2;
        }
        return m;
    }

    void fft(std::vector<std::complex<long double>> &digits, bool invert) {
        for (size_t i = 1, j = 0; i < digits.size(); ++i) {
            size_t bit = digits.size() >> 1;
            for (; j >= bit; bit >>= 1) {
                j -= bit;
            }
            j += bit;
            if (i < j) {
                swap(digits[i], digits[j]);
            }
        }

        for (size_t len = 2; len <= digits.size(); len *= 2) {
            long double ang = 2 * pi / static_cast<long double>(len);
            if (invert)
                ang = -ang;
            std::complex<long double> w_len(cosl(ang), sinl(ang));
            for (size_t i = 0; i < digits.size(); i += len) {
                std::complex<long double> w(1);
                for (size_t j = 0; j < len / 2; ++j) {
                    std::complex<long double> u = digits[i + j];
                    std::complex<long double> v = digits[i + j + len / 2] * w;
                    digits[i + j] = u + v;
                    digits[i + j + len / 2] = u - v;
                    w *= w_len;
                }
            }
        }
        if (invert) {
            long double n = digits.size();
            for (size_t i = 0; i < digits.size(); ++i) {
                digits[i] /= n;
            }
        }
    }
}

class BigInteger;
BigInteger multiply(BigInteger, long long);

class BigInteger {
private:

    static const long long REAL_BASE = 1e3;
    static const long long USER_BASE = 10;
    static const size_t DIGIT_SIZE = 3;

    std::vector<long long> digits;
    bool is_negative;

    void deleteZeroes() {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
        if (!*this) {
            is_negative = false;
        }
    }

    void toCarry() {
        long long carry = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            digits[i] += carry;
            carry = digits[i] / REAL_BASE;
            digits[i] %= REAL_BASE;
        }
        if (carry) {
            digits.push_back(carry);
        }
    }

public:

    void multiply(long long x) {
        is_negative = false;
        for (size_t i = 0; i < digits.size(); ++i) {
            digits[i] *= x;
        }
        toCarry();
        deleteZeroes();
    }

    BigInteger(): digits(1), is_negative(false) {}
    BigInteger(long long x) : is_negative(x < 0) {
        if (x == 0) {
            digits.resize(1);
            return;
        }
        if (is_negative) {
            x = -x;
        }
        while (x) {
            digits.push_back(x % REAL_BASE);
            x /= REAL_BASE;
        }
    }

    void swap(BigInteger& x) {
        std::swap(digits, x.digits);
        std::swap(is_negative, x.is_negative);
    }

    explicit operator bool() const {
        return digits.size() != 1 || digits[0] != 0;
    }

    void applyAbs() {
        is_negative = false;
    }
    void changeSign() {
        if (*this) {
            is_negative = !is_negative;
        }
    }


    BigInteger operator+() const {
        return *this;
    }
    BigInteger operator-() const {
        BigInteger copy = *this;
        copy.changeSign();
        return copy;
    }

    BigInteger& operator+=(const BigInteger& x) {
        if (x.is_negative == is_negative) {
            if (x.digits.size() > digits.size()) {
                digits.resize(x.digits.size());
            }
            long long carry = 0;
            for (size_t i = 0; i < x.digits.size(); ++i) {
                digits[i] += x.digits[i] + carry;
                carry = digits[i] >= REAL_BASE;
                if (carry) {
                    digits[i] -= REAL_BASE;
                }
            }

            if (carry) {
                size_t i = x.digits.size();
                while (i < digits.size() && digits[i] + 1 == REAL_BASE) {
                    digits[i] = 0;
                    ++i;
                }
                if (i == digits.size()) {
                    digits.push_back(0);
                }
                ++digits[i];
            }
        } else {
            if (x.digits.size() >= digits.size()) {
                digits.resize(x.digits.size());
                size_t i = x.digits.size() - 1;
                for (; i > 0 && digits[i] == x.digits[i]; --i) {
                    digits.pop_back();
                }
                is_negative ^= digits[i] < x.digits[i];
            }

            long long carry = 0;
            long long sign = 1 - 2 * (is_negative == x.is_negative);
            for (size_t i = 0; i < x.digits.size(); ++i) {
                digits[i] = sign * (digits[i] - x.digits[i]) - carry;
                carry = digits[i] < 0;
                if (carry) {
                    digits[i] += REAL_BASE;
                }
            }
            if (carry) {
                size_t i = x.digits.size();
                while (digits[i] == 0) {
                    digits[i] = REAL_BASE - 1;
                    ++i;
                }
                --digits[i];
            }
            deleteZeroes();
        }
        return *this;
    }

    BigInteger& operator++() {
        return *this += 1;
    }
    BigInteger operator++(int) {
        BigInteger tmp = *this;
        ++*this;
        return tmp;
    }

    BigInteger& operator-=(const BigInteger& x) {
        changeSign();
        *this += x;
        changeSign();
        return *this;
    }
    BigInteger& operator--() {
        return *this += -1;
    }
    BigInteger operator--(int) {
        BigInteger tmp = *this;
        --*this;
        return tmp;
    }

    BigInteger& operator*=(const BigInteger& x) {
        size_t n = to_pow2(digits.size() + x.digits.size());
        std::vector<std::complex<long double>> fa(n);
        std::vector<std::complex<long double>> fb(n);
        for (size_t i = 0; i < digits.size(); ++i) {
            fa[i] = digits[i];
        }
        for (size_t i = 0; i < x.digits.size(); ++i) {
            fb[i] = x.digits[i];
        }

        fft(fa, false);
        fft (fb, false);
        for (size_t i = 0; i < n; ++i) {
            fa[i] *= fb[i];
        }
        fft (fa, true);

        digits.resize(n);
        for (size_t i = 0; i < n; ++i) {
            digits[i] = static_cast<long long>(roundl(fa[i].real()));
        }

        toCarry();
        is_negative ^= x.is_negative;
        deleteZeroes();
        return *this;
    }

    void multiplyPow10(size_t q) {
        if (*this) {
            digits.insert(digits.begin(), q / DIGIT_SIZE, 0);
            if (q %= DIGIT_SIZE) {
                long long x = USER_BASE;
                while (--q) {
                    x *= USER_BASE;
                }
                multiply(x);
            }
        }
    }

    std::pair<BigInteger, BigInteger> div_mod(const BigInteger& x) const {
        BigInteger div;
        BigInteger mod;
        for (size_t i = digits.size() - 1; i < digits.size(); --i) {
            if (mod.digits.back()) {
                mod.digits.insert(mod.digits.begin(), 1, digits[i]);
            } else {
                mod.digits[0] = digits[i];
            }
            div.digits.push_back(0);

            if (mod >= x) {
                long long &left = div.digits.back();
                long long right = REAL_BASE;
                while (left + 1 < right) {
                    long long md = (left + right) / 2;
                    if (::multiply(x, md) <= mod) {
                        left = md;
                    } else {
                        right = md;
                    }
                }
                mod -= ::multiply(x, left);
            }
        }
        reverse(div.digits.begin(), div.digits.end());
        div.is_negative = is_negative ^ x.is_negative;
        mod.is_negative = mod && is_negative;
        div.deleteZeroes();
        return {div, mod};
    }

    BigInteger& operator/=(const BigInteger& x) {
        return *this = div_mod(x).first;
    }

    BigInteger& operator%=(const BigInteger& x) {
        return *this = div_mod(x).second;
    }


    bool operator==(const BigInteger&) const = default;
    friend std::strong_ordering operator<=>(const BigInteger& x, const BigInteger& y) {
        std::strong_ordering less = std::strong_ordering::less;
        std::strong_ordering equal = std::strong_ordering::equal;
        std::strong_ordering greater = std::strong_ordering::greater;

        if (x.is_negative && !y.is_negative) {
            return less;
        }
        if (!x.is_negative && y.is_negative) {
            return greater;
        }
        if (x.is_negative && y.is_negative) {
            std::swap(less, greater);
        }
        if (x.digits.size() < y.digits.size()) {
            return less;
        }
        if (x.digits.size() > y.digits.size()) {
            return greater;
        }

        size_t i = x.digits.size() - 1;
        while (i < x.digits.size() && x.digits[i] == y.digits[i]) {
            --i;
        }
        if (i > x.digits.size()) {
            return equal;
        }

        return x.digits[i] < y.digits[i] ? less : greater;
    }

    std::string toString() const {
        std::string s(DIGIT_SIZE * digits.size(), '\0');
        for (size_t i = 0; i < digits.size(); ++i) {
            long long t = digits[i];
            for (size_t j = DIGIT_SIZE * i; j < DIGIT_SIZE * (i + 1); ++j) {
                s[j] = std::to_string(t % USER_BASE)[0];
                t /= USER_BASE;
            }
        }
        while (s.size() > 1 && s.back() == '0') {
            s.pop_back();
        }
        if (is_negative) {
            s.push_back('-');
        }
        reverse(s.begin(), s.end());
        return s;
    }

    friend std::istream& operator>>(std::istream& in, BigInteger& x) {
        std::string s;
        in >> s;
        reverse(s.begin(), s.end());
        if (s.back() == '-') {
            x.is_negative = true;
            s.pop_back();
        } else {
            x.is_negative = false;
        }
        while (s.size() % DIGIT_SIZE) {
            s.push_back('0');
        }
        x.digits.resize(s.size() / DIGIT_SIZE);

        for (size_t i = 0; i < x.digits.size(); ++i) {
            x.digits[i] = 0;
            for (size_t j = DIGIT_SIZE * (i + 1); j > DIGIT_SIZE * i; --j) {
                x.digits[i] = USER_BASE * x.digits[i] + (s[j - 1] - '0');
            }
        }
        x.deleteZeroes();
        return in;
    }

};

BigInteger multiply(BigInteger ans, long long x) {
    ans.multiply(x);
    return ans;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& x) {
    return out << x.toString();
}
BigInteger operator"" _bi(unsigned long long n) {
    return BigInteger(static_cast<long long>(n));
}

BigInteger operator+(BigInteger x, const BigInteger& y) {
    x += y;
    return x;
}
BigInteger operator-(BigInteger x, const BigInteger& y) {
    x -= y;
    return x;
}
BigInteger operator*(BigInteger x, const BigInteger& y) {
    x *= y;
    return x;
}
BigInteger operator/(BigInteger x, const BigInteger& y) {
    x /= y;
    return x;
}
BigInteger operator%(BigInteger x, const BigInteger& y) {
    x %= y;
    return x;
}


BigInteger gcd(BigInteger a, BigInteger b) {
    a.applyAbs();
    b.applyAbs();

    while (b) {
        a %= b;
        a.swap(b);
    }
    return a;
}

class Rational {
private:
    static const size_t MANTISSA_SIZE = 20;
    BigInteger x;
    BigInteger y;

    void reduce() {
        if (y < 0) {
            x.changeSign();
            y.changeSign();
        }
        BigInteger g = gcd(x, y);
        x /= g;
        y /= g;
    }
public:
    Rational(const BigInteger& x): x(x), y(1) {}
    Rational(const BigInteger& x, const BigInteger& y): x(x), y(y) {
        reduce();
    }

    Rational(long long x) : x(x), y(1) {}
    Rational(long long x, long long y) : x(x), y(y) {
        reduce();
    }
    Rational() : x(0), y(1) {}

    explicit operator double() const {
        return atof(asDecimal(MANTISSA_SIZE).c_str());
    }

    Rational operator+() const {
        return *this;
    }
    Rational operator-() const {
        Rational copy(*this);
        copy.x.changeSign();
        return copy;
    }

    Rational& operator+=(const Rational& t) {
        x *= t.y;
        x += t.x * y;
        y *= t.y;
        reduce();
        return *this;
    }
    Rational& operator-=(const Rational& t) {
        x.changeSign();
        *this += t;
        x.changeSign();
        return *this;
    }

    Rational& operator*=(const Rational& t) {
        x *= t.x;
        y *= t.y;
        reduce();
        return *this;
    }
    Rational& operator/=(const Rational& t) {
        x *= t.y;
        y *= t.x;
        reduce();
        return *this;
    }

    bool operator==(const Rational&) const = default;
    friend std::strong_ordering operator<=>(const Rational& digits, const Rational& b) {
        return digits.x * b.y <=> b.x * digits.y;
    }

    std::string toString() const {
        return x.toString() + (y != 1_bi ? '/' + y.toString() : "");
    }

    std::string asDecimal(size_t precision=0) const {
        auto [div, mod] = x.div_mod(y);
        std::string integer = div.toString();
        if (precision == 0) {
            return integer;
        }
        mod.applyAbs();
        mod.multiplyPow10(precision);
        mod /= y;
        std::string s = mod.toString();
        if (integer == "0" && x < 0 && mod) {
            integer = "-0";
        }
        return integer + '.' + std::string(precision - s.size(), '0') + s;
    }
};

Rational operator+(Rational a, const Rational& b) {
    a += b;
    return a;
}
Rational operator-(Rational a, const Rational& b) {
    a -= b;
    return a;
}
Rational operator*(Rational a, const Rational& b) {
    a *= b;
    return a;
}
Rational operator/(Rational a, const Rational& b) {
    a /= b;
    return a;
}

std::istream& operator>>(std::istream& in, Rational& x) {
    BigInteger t;
    in >> t;
    x = Rational(t);
    return in;
}
std::ostream& operator<<(std::ostream& out, const Rational& x) {
    return out << x.toString();
}

template <size_t i, size_t N>
struct isPrime {
    static const bool value = N % i && isPrime<(i * i <= N) * (i + 1), N>::value;
};

template <size_t N>
struct isPrime<0, N> {
    static const bool value = true;
};
template <>
struct isPrime<2, 2> {
    static const bool value = true;
};

template <size_t N>
class Residue {
public:
    size_t x = 0;
    Residue() = default;

    Residue(int x) :
            x(static_cast<size_t>((x % static_cast<long long>(N) + static_cast<long long>(N))) % N) {}

    operator int() const {
        return x;
    }

    Residue& operator*=(const Residue& t) {
        x = (x * t.x) % N;
        return *this;
    }
    Residue& operator+=(const Residue& t) {
        x += t.x;
        if (x >= N) {
            x -= N;
        }
        return *this;
    }
    Residue& operator-=(const Residue& t) {
        x -= t.x;
        if (x >= N) {
            x += N;
        }
        return *this;
    }

    Residue pow(size_t b) {
        Residue ans = 1;
        Residue a = *this;
        while (b) {
            if (b % 2) {
                ans = ans * a;
            }
            a *= a;
            b /= 2;
        }
        return ans;
    }

    Residue& operator/=(Residue t) {
        static_assert(isPrime<2, N>::value);
        return *this *= t.pow(N - 2);
    }
};

template <size_t N>
Residue<N> operator*(Residue<N> a, const Residue<N>& b) {
    a *= b;
    return a;
}
template <size_t N>
Residue<N> operator+(Residue<N> a, const Residue<N>& b) {
    a += b;
    return a;
}
template <size_t N>
Residue<N> operator-(Residue<N> a, const Residue<N>& b) {
    a -= b;
    return a;
}

template <size_t N>
Residue<N> operator/(Residue<N> a, const Residue<N>& b) {
    a /= b;
    return a;
}

template <size_t N, size_t M>
struct toPow2 {
    static const size_t value = M >= N ? M : toPow2<N, (M < N) * (2 * M)>::value;
};
template <size_t N>
struct toPow2<N, 0> {
    static const size_t value = 0;
};

template <size_t N, size_t M, typename Field=Rational>
class Matrix {
private:
    std::array<std::array<Field, M>, N> data;

public:
    Matrix() {
        if constexpr (!std::is_class_v<Field>) {
            for (size_t i = 0; i < N; ++i) {
                std::fill(data[i].begin(), data[i].end(), Field(0));
            }
        }
    }

    Matrix(const std::initializer_list<std::initializer_list<Field>>& row) {
        auto t = data.begin();
        for (auto& it : row) {
            std::copy(it.begin(), it.end(), t->begin());
            ++t;
        }
    }

    bool operator==(const Matrix<N, M, Field>& x) const {
        return data == x.data;
    }

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& t) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                data[i][j] += t.data[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& t) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                data[i][j] -= t.data[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field>& operator*=(const Field& k) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                data[i][j] *= k;
            }
        }
        return *this;
    }

    template <size_t K>
    Matrix<N, M, Field>& operator*=(const Matrix<M, K, Field>& x) {
        return *this = *this * x;
    }

    const std::array<Field, M>& operator[](size_t i) const {
        return data[i];
    }

    std::array<Field, M>& operator[](size_t i) {
        return data[i];
    }


    Matrix<M, N, Field> transposed() const {
        Matrix<M, N, Field> ans;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                ans[i][j] = data[j][i];
            }
        }
        return ans;
    }

    std::pair<size_t, bool> toRowEchelonForm() {
        size_t rank = 0;
        bool parity = false;
        for (size_t i = 0, k = 0; i < N && k < M; ++k) {
            if (data[i][k] == Field(0)) {
                size_t j = i + 1;
                while (j < N && data[j][k] == Field(0)) {
                    ++j;
                }
                if (j == N) {
                    continue;
                }
                std::swap(data[i], data[j]);
                parity ^= true;
            }

            for (size_t j = i + 1; j < N; ++j) {
                Field t = data[j][k] / data[i][k];
                for (size_t c = k; c < M; ++c) {
                    data[j][c] -= t * data[i][c];
                }
            }
            ++rank;
            ++i;
        }
        return {rank, parity};
    }

    Field det() const {
        static_assert(N == M);
        Matrix copy(*this);
        auto [rank, parity] = copy.toRowEchelonForm();
        if (rank < N) {
            return Field(0);
        }
        Field ans = parity ? -1 : 1;
        for (size_t i = 0; i < N; ++i) {
            ans *= copy.data[i][i];
        }
        return ans;
    }


    size_t rank() const {
        Matrix copy(*this);
        return copy.toRowEchelonForm().first;
    }

    void invert() {
        static_assert(N == M);
        Matrix ans;
        for (size_t i = 0; i < N; ++i) {
            ans[i][i] = Field(1);
        }
        for (size_t i = 0; i < N; ++i) {
            size_t j = i;
            while (j < N && data[j][i] == Field(0)) {
                ++j;
            }
            if (j == N) {
                throw "Singular matrix inversion";
            }
            std::swap(data[i], data[j]);
            std::swap(ans.data[i], ans.data[j]);

            for (j = 0; j < N; ++j) {
                ans.data[i][j] /= data[i][i];
            }
            for (j = i + 1; j < N; ++j) {
                data[i][j] /= data[i][i];
            }
            data[i][i] = Field(1);

            for (j = i + 1; j < N; ++j) {
                for (size_t k = i + 1; k < N; ++k) {
                    data[j][k] -= data[j][i] * data[i][k];
                }
                for (size_t k = 0; k < N; ++k) {
                    ans.data[j][k] -= data[j][i] * ans.data[i][k];
                }
                data[j][i] = Field(0);
            }
        }

        for (size_t i = N - 1; i < N; --i) {
            for (size_t j = 0; j < i; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    ans.data[j][k] -= data[j][i] * ans.data[i][k];
                }
            }
        }
        std::swap(data, ans.data);
    }

    Matrix<N, N, Field> inverted() const {
        Matrix ans(*this);
        ans.invert();
        return ans;
    }

    Field trace() const {
        size_t n = std::min(N, M);
        Field ans = 0;
        for (size_t i = 0; i < n; ++i) {
            ans += data[i][i];
        }
        return ans;
    }

    std::array<Field, M> getRow(size_t i) {
        return data[i];
    }

    std::array<Field, N> getColumn(size_t i) {
        std::array<Field, N> ans;
        for (size_t j = 0; j < N; ++j) {
            ans[j] = data[j][i];
        }
        return ans;
    }
};

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator+(Matrix<N, M, Field> a, const Matrix<N, M, Field>& b) {
    a += b;
    return a;
}
template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator-(Matrix<N, M, Field> a, const Matrix<N, M, Field>& b) {
    a -= b;
    return a;
}
template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(Matrix<N, M, Field> a, const Field& b) {
    a *= b;
    return a;
}
template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Field& a, Matrix<N, M, Field> b) {
    b *= a;
    return b;
}

namespace {
    template<size_t N, typename Field>
    Matrix<N, N, Field> multiply(Matrix<N, N, Field>& a, Matrix<N, N, Field>& b) {
        Matrix<N, N, Field> ans;
        if (N == 1) {
            ans[0][0] = a[0][0] * b[0][0];
            return ans;
        }
        const size_t M = N / 2;
        Matrix<M, M, Field> a11;
        Matrix<M, M, Field> a12;
        Matrix<M, M, Field> a21;
        Matrix<M, M, Field> a22;

        Matrix<M, M, Field> b11;
        Matrix<M, M, Field> b12;
        Matrix<M, M, Field> b21;
        Matrix<M, M, Field> b22;

        for (size_t i = 0; i < M; ++i) {
            std::copy(a[i].begin(), a[i].begin() + M, a11[i].begin());
            std::copy(a[i].begin() + M, a[i].end(), a12[i].begin());
            std::copy(a[i + M].begin(), a[i + M].begin() + M, a21[i].begin());
            std::copy(a[i + M].begin() + M, a[i + M].end(), a22[i].begin());

            std::copy(b[i].begin(), b[i].begin() + M, b11[i].begin());
            std::copy(b[i].begin() + M, b[i].end(), b12[i].begin());
            std::copy(b[i + M].begin(), b[i + M].begin() + M, b21[i].begin());
            std::copy(b[i + M].begin() + M, b[i + M].end(), b22[i].begin());
        }

        Matrix<M, M, Field> s1 = a21 + a22;
        Matrix<M, M, Field> s2 = s1 - a11;
        Matrix<M, M, Field> s3 = a11 - a21;
        Matrix<M, M, Field> s4 = a12 - s2;
        Matrix<M, M, Field> s5 = b12 - b11;
        Matrix<M, M, Field> s6 = b22 - s5;
        Matrix<M, M, Field> s7 = b22 - b12;
        Matrix<M, M, Field> s8 = s6 - b21;

        Matrix<M, M, Field> p1 = multiply(s2, s6);
        Matrix<M, M, Field> p2 = multiply(a11, b11);
        Matrix<M, M, Field> p3 = multiply(a12, b21);
        Matrix<M, M, Field> p4 = multiply(s3, s7);
        Matrix<M, M, Field> p5 = multiply(s1, s5);
        Matrix<M, M, Field> p6 = multiply(s4, b22);
        Matrix<M, M, Field> p7 = multiply(a22, s8);

        Matrix<M, M, Field> t1 = p1 + p2;
        Matrix<M, M, Field> t2 = t1 + p4;

        Matrix<M, M, Field> c11 = p2 + p3;
        Matrix<M, M, Field> c12 = t1 + p5 + p6;
        Matrix<M, M, Field> c21 = t2 - p7;
        Matrix<M, M, Field> c22 = t2 + p5;

        for (size_t i = 0; i < M; ++i) {
            std::copy(c11[i].begin(), c11[i].end(), ans[i].begin());
            std::copy(c12[i].begin(), c12[i].end(), ans[i].begin() + M);
            std::copy(c21[i].begin(), c21[i].end(), ans[i + M].begin());
            std::copy(c22[i].begin(), c22[i].end(), ans[i + M].begin() + M);
        }
        return ans;
    }
}

template <size_t N, size_t M, size_t K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field>& x, const Matrix<M, K, Field>& y) {
    Matrix<N, K, Field> ans;
    const size_t n = toPow2<std::max({N, M, K}), 1>::value;
    Matrix<n, n, Field> data;
    Matrix<n, n, Field> b;
    for (size_t i = 0; i < N; ++i) {
        std::copy(x[i].begin(), x[i].end(), data[i].begin());
    }
    for (size_t i = 0; i < M; ++i) {
        std::copy(y[i].begin(), y[i].end(), b[i].begin());
    }
    Matrix<n, n, Field> c = multiply(data, b);
    for (size_t i = 0; i < N; ++i) {
        std::copy(c[i].begin(), c[i].begin() + K, ans[i].begin());
    }
    return ans;
}

template <size_t N, typename T=Rational>
using SquareMatrix = Matrix<N, N, T>;
