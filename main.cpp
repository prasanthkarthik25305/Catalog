#include <bits/stdc++.h>
using namespace std;

// ======= Big Integer string utilities (with sign support) =======
string stripZeros(string s) {
    int i = 0;
    while (i < (int)s.length() - 1 && s[i] == '0') i++;
    return s.substr(i);
}

// Compare absolute (a,b): returns -1, 0, 1
int cmpAbs(const string &a, const string &b) {
    string sa = stripZeros(a), sb = stripZeros(b);
    if (sa.size() != sb.size()) return sa.size() < sb.size() ? -1 : 1;
    if (sa == sb) return 0;
    return sa < sb ? -1 : 1;
}

// Add two non-negative int strings
string addUnsigned(string a, string b) {
    if (a.size() < b.size()) swap(a, b);
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    string res;
    int carry = 0;
    for (size_t i = 0; i < a.size(); i++) {
        int sum = (a[i] - '0') + carry + (i < b.size() ? (b[i] - '0') : 0);
        carry = sum / 10;
        res.push_back((sum % 10) + '0');
    }
    if (carry) res.push_back(carry + '0');
    reverse(res.begin(), res.end());
    return stripZeros(res);
}

// Subtract (a >= b), unsigned
string subUnsigned(string a, string b) {
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    string res;
    int carry = 0;
    for (size_t i = 0; i < a.size(); i++) {
        int diff = (a[i] - '0') - carry - (i < b.size() ? (b[i] - '0') : 0);
        if (diff < 0) { diff += 10; carry = 1; } else carry = 0;
        res.push_back(diff + '0');
    }
    while (res.size() > 1 && res.back() == '0') res.pop_back();
    reverse(res.begin(), res.end());
    return stripZeros(res);
}

// Addition for signed big int strings
string addBig(string a, string b) {
    bool negA = false, negB = false;
    if (a[0] == '-') { negA = true; a = a.substr(1); }
    if (b[0] == '-') { negB = true; b = b.substr(1); }
    if (negA == negB) {
        string sum = addUnsigned(a, b);
        return (negA && sum != "0") ? "-" + sum : sum;
    } else {
        int cmp = cmpAbs(a, b);
        if (cmp == 0) return "0";
        if (cmp > 0) { // a > b
            string diff = subUnsigned(a, b);
            return (negA ? "-" + diff : diff);
        } else {
            string diff = subUnsigned(b, a);
            return (negB ? "-" + diff : diff);
        }
    }
}

// Subtraction for signed big int strings: a - b
string subBig(string a, string b) {
    bool negB = (b[0] == '-');
    b = negB ? b.substr(1) : "-" + b;
    return addBig(a, b);
}

// Multiply unsigned
string mulUnsigned(string a, string b) {
    vector<int> res(a.size() + b.size(), 0);
    for (int i = a.size() - 1; i >= 0; i--) {
        for (int j = b.size() - 1; j >= 0; j--) {
            int mul = (a[i] - '0') * (b[j] - '0');
            int sum = mul + res[i + j + 1];
            res[i + j + 1] = sum % 10;
            res[i + j] += sum / 10;
        }
    }
    string s;
    for (int num : res) {
        if (!(s.empty() && num == 0)) s.push_back(num + '0');
    }
    return s.empty() ? "0" : s;
}

// Signed multiply
string mulBig(string a, string b) {
    bool neg = false;
    if (a[0] == '-') { neg = !neg; a = a.substr(1); }
    if (b[0] == '-') { neg = !neg; b = b.substr(1); }
    string prod = mulUnsigned(a, b);
    if (neg && prod != "0") prod = "-" + prod;
    return prod;
}

// Integer division (b fits in long long), signed
string divBig(string a, string b) {
    bool neg = false;
    if (a[0] == '-') { neg = !neg; a = a.substr(1); }
    if (b[0] == '-') { neg = !neg; b = b.substr(1); }
    long long divisor = stoll(b);
    string res;
    long long cur = 0;
    for (char c : a) {
        cur = cur * 10 + (c - '0');
        if (!res.empty() || cur / divisor > 0) {
            res.push_back((cur / divisor) + '0');
        }
        cur %= divisor;
    }
    if (res.empty()) res = "0";
    return (neg && res != "0") ? "-" + res : res;
}

// Convert from any base to decimal string
string baseToDec(const string &val, int base) {
    string result = "0";
    string baseStr = to_string(base);
    for (char c : val) {
        int digit;
        if (isdigit(c)) digit = c - '0';
        else digit = tolower(c) - 'a' + 10;
        result = mulBig(result, baseStr);
        result = addBig(result, to_string(digit));
    }
    return result;
}

bool isZero(const string &s) {
    return s == "0" || s == "-0";
}

// ---------- Gaussian elimination for big integer strings ----------
vector<string> gaussianSolve(vector<vector<string>> mat) {
    int n = mat.size();
    int m = mat[0].size() - 1; // number of unknowns = k+1

    int row = 0;
    for (int col = 0; col < m && row < n; col++) {
        // Find pivot
        int pivot = row;
        for (int i = row; i < n; i++) {
            if (!isZero(mat[i][col])) { pivot = i; break; }
        }
        if (isZero(mat[pivot][col])) continue;

        if (pivot != row) swap(mat[pivot], mat[row]);

        // Normalize pivot row (divide entire row by mat[row][col])
        string pivotVal = mat[row][col];
        for (int j = col; j <= m; j++)
            mat[row][j] = divBig(mat[row][j], pivotVal); // trust integer division exactness

        // Eliminate in other rows
        for (int i = 0; i < n; i++) {
            if (i != row && !isZero(mat[i][col])) {
                string factor = mat[i][col];
                for (int j = col; j <= m; j++)
                    mat[i][j] = subBig(mat[i][j], mulBig(factor, mat[row][j]));
            }
        }
        row++;
    }

    // Extract solution (last column after elimination)
    vector<string> sol(m, "0");
    for (int i = 0; i < m; i++) sol[i] = mat[i][m];
    return sol;
}

// ----------- Main program -----------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string json, line;
    while (getline(cin, line)) json += line;

    // Parse "n" and "k"
    int n = 0, k = 0;
    {
        regex nk("\"n\"\\s*:\\s*(\\d+).*?\"k\"\\s*:\\s*(\\d+)");
        smatch m;
        if (regex_search(json, m, nk)) {
            n = stoi(m[1]);
            k = stoi(m[2]);
        }
    }

    // Parse all points: key: { "base": "...", "value": "..." }
    map<int, pair<int, string>> data;
    regex point("\"(\\d+)\"\\s*:\\s*\\{\\s*\"base\"\\s*:\\s*\"(\\d+)\"\\s*,\\s*\"value\"\\s*:\\s*\"([^\"]+)\"");
    smatch m;
    string::const_iterator start(json.cbegin());
    while (regex_search(start, json.cend(), m, point)) {
        int x = stoi(m[1]);
        int b = stoi(m[2]);
        string val = m[3];
        data[x] = {b, val};
        start = m.suffix().first;
    }

    // Build augmented matrix for degree k polynomial a_k*x^k + ... + a_0 = y
    int sz = k + 1;
    vector<vector<string>> mat(sz, vector<string>(sz + 1, "0"));

    int idx = 0;
    for (auto &p : data) {
        if (idx >= sz) break; // use first k+1 points for interpolation
        string xi = to_string(p.first);
        string yi = baseToDec(p.second.second, p.second.first);

        // Fill powers of x
        string powx = "1";
        for (int j = sz - 1; j >= 0; j--) { // x^j
            mat[idx][sz - 1 - j] = powx;
            powx = mulBig(powx, xi);
        }
        mat[idx][sz] = yi;
        idx++;
    }

    // Solve
    vector<string> coeff = gaussianSolve(mat);

    // Output constant term (a0 at the end)
    // coeff[0]: x^k, ..., coeff[k]: x^0 constant term
    cout << coeff[k] << "\n";

    return 0;
}

