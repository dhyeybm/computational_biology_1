#include <bits/stdc++.h>
#include <chrono>
#include <cstdlib>
#include <random>
#include <iostream>

using namespace std;
using namespace std::chrono;


typedef long long ll;

typedef vector<int> vint;
typedef pair<vint, vint> vpair;
typedef pair<string, string> pstring;

ulong vector_to_index(vector<char>::iterator it, vector<char>::iterator end, int base, int offset) {
    ulong result = 0;
    while (it != end) {
        result = result * base + ((*it) - offset);
        it++;
    }
    return result;
}

ulong vector_to_index(vector<int>::iterator it, vector<int>::iterator end, int base, int offset) {
    ulong result = 0;
    while (it != end) {
        result = result * base + ((*it) - offset);
        it++;
    }
    return result;
}

inline ulong vector_to_index(const char *arr, int n, int base, int offset) {
    ulong result = 0, i = 0;
    while (i < n) {
        result = result * base + (arr[i] - offset);
        i += 1;
    }
    return result;
}
inline ulong vector_to_index(const int *arr, int n, int base, int offset) {
    ulong result = 0, i = 0;
    while (i < n) {
        result = result * base + (arr[i] - offset);
        i += 1;
    }
    return result;
}

ulong int_pow(int x, int m) {
    ulong result = 1;
    for (int i = 0; i < m; i++) {
        result *= x;
    }
    return result;
}

char *load_cache(const string &filename, int t, int s) {
    ifstream f(filename, ios::binary | ios::ate);
    streamsize size = f.tellg();
    f.seekg(0, ios::beg);

    ulong max_value = int_pow(s, t);
    ulong max_value2 = int_pow(3, t);
    ulong total_size = 2 * max_value * max_value * max_value2 * max_value2 * t;

    if (total_size != size) {
        cerr << "Cache size does not match expected size" << endl;
        exit(-1);
    }

    cerr << "Array size: " << size << endl;

    char *cache = new char[size];
    if (f.read(cache, size)) {
        return cache;
    } else {
        fprintf(stderr, "Failed to read file");
        exit(1);
    }
}

void fill_partial_edit_distance(const vector<char> &b, const vector<char> &c,
                                const vector<char> &v, const vector<char> &w,
                                vector<char> &d, vector<char> &e) {
    /*
     * Some other implicit restrictions
     * b.size() == v.size()
     * c.size() == w.size()
     * d.size() == b.size()
     * e.size() == c.size()
     */
    unsigned long m = v.size() + 1, n = w.size() + 1;
    vector<vector<char> > dp(m, vector<char>(n));

    // Initializations
    dp[0][0] = 0;

    // Initialize left border of dp matrix
    for (int i = 1; i < m; i++) {
        dp[i][0] = dp[i - 1][0] + b[i - 1];
    }

    // Initialize top-border of dp matrix
    for (int j = 1; j < n; j++) {
        dp[0][j] = dp[0][j - 1] + c[j - 1];
    }

    // Fill interior of dp tale
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            if (v[i - 1] == w[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min(dp[i - 1][j], min(dp[i][j - 1], dp[i - 1][j - 1])) + 1;
            }
        }
    }

    for (int i = 0; i < d.size(); i++) {
        d[i] = dp[i + 1][n - 1];
    }

    for (int j = 0; j < e.size(); j++) {
        e[j] = dp[m - 1][j + 1];
    }
}

void fill_base_values(vector<char> &v, ulong x, int base, int offset) {
    int i = v.size() - 1;
    while (i >= 0) {
        v[i] = char(x % base) + offset;
        x /= base;
        i -= 1;
    }
}

char *cache_values(int t, int s) {
    ulong tpow = int_pow(s, t);
    ulong spow = int_pow(3, t);

    ulong max_value = int_pow(s, t);
    ulong max_value2 = int_pow(3, t);

    ulong total_size = 2 * tpow * spow * max_value * max_value2 * t;
    char *cache = new char[total_size];

    cerr << "Total size: " << total_size << endl;
    char *cache_iter = cache;
    vector<char> v(t), w(t), b(t), c(t), d(t), e(t);

    for (ulong x = 0; x < max_value; x++) {
        cerr << "Completed |   x: " << x + 1 << "/" << max_value << endl;
        fill_base_values(v, x, s, 0);
        for (ulong y = 0; y < max_value; y++) {
            fill_base_values(w, y, s, 0);

            for (ulong bx = 0; bx < max_value2; bx++) {
                fill_base_values(b, bx, 3, -1);

                for (ulong cx = 0; cx < max_value2; cx++) {
                    fill_base_values(c, cx, 3, -1);

                    fill_partial_edit_distance(b, c, v, w, d, e);

                    for (int i = 0; i < t; i++) {
                        *cache_iter = d[i];
                        cache_iter++;
                    }
                    for (int i = 0; i < t; i++) {
                        *cache_iter = e[i];
                        cache_iter++;
                    }
                }
            }
        }
    }

    return cache;
}

inline bool file_exists(const string &name) {
    ifstream f(name.c_str());
    return f.good();
}

char *calculate_or_load_cache(int t, int s) {
    string filename = "_cache_" + to_string(t) + "_" + to_string(s);
    if (!file_exists(filename)) {
        fprintf(stderr, "Creating cache for t=%d, s=%d\n", t, s);
        char *cache = cache_values(t, s);
        ulong max_value = int_pow(s, t);
        ulong max_value2 = int_pow(3, t);
        ulong total_size = 2 * max_value * max_value * max_value2 * max_value2 * t;
        ofstream f(filename, ios::binary);
        f.write(cache, total_size);
        f.close();
        return cache;
    } else {
        fprintf(stderr, "Loading cache from file %s\n", filename.c_str());
        return load_cache(filename, t, s);
    }
}


class Encoder {
private:
    map<char, char> e;
    map<char, char> d;
public:
    Encoder(string s) {
        char present[256];
        for (int i = 0; i < 256; i++) { present[i] = false; }
        for (int i = 0; i < s.length(); i++) { present[s[i]] = true; }

        int count = 0;
        for (int i = 0; i < 256; i++) {
            if (present[i]) {
                e[(char) i] = (char) ('0' + count);
                d[(char) ('0' + count)] = (char) i;
                count++;
            }
        }
    }

    size_t charset_size() {
        return e.size();
    }

    string encode(string s) {
        vector<char> v;
        for (int i = 0; i < s.length(); i++) {
            v.push_back(e[s[i]]);
        }
        return string(v.begin(), v.end());
    }

    string decode(string s) {
        vector<char> v;
        for (int i = 0; i < s.length(); i++) {
            v.push_back(s[i] == '-' ? '-' : d[s[i]]);
        }
        return string(v.begin(), v.end());
    }
};

vector<int> edit_distance_rightmost_col(const vector<char> &v, const vector<char> &w) {
    unsigned long m = v.size(), n = w.size();
    vector<vector<int> > dp(m, vector<int>(n));
    vector<int> result(m);

    // Initialize borders
    for (int i = 0; i < m; i++) {
        dp[i][0] = i;
    }

    for (int i = 0; i < n; i++) {
        dp[0][i] = i;
    }

    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            if (v[i] == w[j]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min(dp[i - 1][j], min(dp[i][j - 1], dp[i - 1][j - 1])) + 1;
            }
        }
    }

    for (int i = 0; i < m; i++) {
        result[i] = dp[i][n - 1];
    }
    return result;
}

vector<int> _russians(const vector<char> &v, const vector<char> &w, int s, int t, const char *cache) {
    ulong tpow = int_pow(s, t);
    ulong spow = int_pow(3, t);

    ulong m = v.size();
    ulong n = w.size();

    char b[6];
    char c[6];

    // Initialize
    for (int i = 0; i < t; i++) {
        b[i] = 1;
        c[i] = 1;
    }

    vector<int> rightmost_col(v.size());
    vector<int> row_above(w.size()), row_below(w.size());

    for (int j = 0; j < n; j++) {
        row_above[j] = j;
    }
    rightmost_col[0] = w.size() - 1;

    vector<int> w_index(w.size() / t);
    for (int j=0; j + t < n; j += t) {
        w_index[j / t] = (int)vector_to_index(&w[j + 1], t, s, 0);
    }

    for (int i = 0; i + t < m; i += t) {
        // Starting new row. Moving down a block. Update b
        for (int k = 0; k < t; k++) {
            b[k] = 1;
        }
        ulong base_index = vector_to_index(&v[i + 1], t, s, 0);
        ulong index;
        for (int j = 0; j + t < n; j += t) {
            // Moving right. Update c
            for (int k = 0; k < t; k++) {
                c[k] = (char)(row_above[j + k + 1] - row_above[j + k]);
            }

            // Get index into cache
            index = base_index;
//            index = index * tpow + vector_to_index(&w[j + 1], t, s, 0);
            index = index * tpow + w_index[j / t];
            index = index * spow + vector_to_index(&b[0], t, 3, -1);
            index = index * spow + vector_to_index(&c[0], t, 3, -1);
            index = index * 2 * t;

            for (int k = 0; k < t; k++) {
                if (k == 0) {
                    b[0] = row_above[j] + cache[index] - row_above[j + t];
                } else {
                    b[k] = cache[index + k] - cache[index + k - 1];
                }
                if (j + 2 * t >= n) {
                    rightmost_col[i + k + 1] = row_above[n - 1 - t] + cache[index + k];
                }
            }

            index += t;
            for (int k = 0; k < t; k++) {
                row_below[j + k + 1] = row_above[j] + cache[index + k];
            }
        }

        // Starting new row after this. So row below will be the row above
        row_above.swap(row_below);
        row_above[0] = i + t;
    }
    return rightmost_col;
}


vector<int> russians(vector<char> v, vector<char> w, int s, int t, const char *cache) {
    if (v.size() >= t + 1 && w.size() >= t + 1) {
        // We hack our way through to get russians to work for different sizes
        ulong initial_v_size = v.size();
        ulong w2_size;
        if (w.size() % t == 0) {
            w2_size = w.size() - t + 1;
        } else {
            w2_size = w.size() - w.size() % t + 1;
        }
        while (v.size() % t != 1) {
            v.push_back(0);
        }

        vector<char> w2(w.begin(), w.begin() + w2_size);
        vector<int> result_temp = _russians(v, w2, s, t, cache);
        vector<int> result = result_temp;

        // Extend this edit distance
        for (int j = w2.size(); j < w.size(); j++) {
            result[0] = j;
            for (int i = 1; i < v.size(); i++) {
                if (v[i] == w[j]) {
                    result[i] = result_temp[i - 1];
                } else {
                    result[i] = min(result_temp[i], min(result_temp[i - 1], result[i - 1])) + 1;
                }
            }
            result_temp = result;
        }
        return vector<int>(result.begin(), result.begin() + initial_v_size);
    } else {
        return edit_distance_rightmost_col(v, w);
    }
}

vector<int> russians_prefix(int i1, int j1, int i2, int j2,
                            string &seq1, string &seq2, char *cache, int t, int s) {
    vector<char> v, w;
    for (int i = i1; i <= i2; i++) {
        v.push_back(seq1[i] - '0');
    }
    for (int j = j1; j <= j2; j++) {
        w.push_back(seq2[j] - '0');
    }

    vector<int> result = russians(v, w, s, t, cache);
    for (int i = 0; i < result.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}

vector<int> russians_suffix(int i1, int j1, int i2, int j2,
                            string &seq1, string &seq2, char *cache, int t, int s) {
    vector<char> v, w;
    for (int i = i1 + 1; i <= i2; i++) {
        v.push_back(seq1[i] - '0');
    }
    v.push_back(0);
    reverse(v.begin(), v.end());

    for (int j = j1 + 1; j <= j2; j++) {
        w.push_back(seq2[j] - '0');
    }
    w.push_back(0);
    reverse(w.begin(), w.end());

    vector<int> result = russians(v, w, s, t, cache);
    reverse(result.begin(), result.end());
    for (int i = 0; i < result.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}

vector<int> _find_i_star_list_russian(int i1, int j1, int i2, int j2, int j, string &seq1, string &seq2,
                                      char *cache, int t, int s) {
    vector<int> pre = russians_prefix(i1, j1, i2, j, seq1, seq2, cache, t, s);
    vector<int> suf = russians_suffix(i1, j, i2, j2, seq1, seq2, cache, t, s);

    int max = pre[0] + suf[0], list_begin_ix = 0;
    for (int i = 0; i < pre.size(); ++i) {
        if (max < pre[i] + suf[i]) {
            max = pre[i] + suf[i];
            list_begin_ix = i;
        }
    }
    vector<int> i_star_list;
    i_star_list.push_back(i1 + list_begin_ix);
    for (int i = list_begin_ix + 1; i < pre.size(); ++i) {
        if ((pre[i] == pre[i - 1] - 1) && max == pre[i] + suf[i])
            i_star_list.push_back(i1 + i);
        else
            break;
    }

    return i_star_list;
}

void _hirschberg_russians(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring,
                          map<int, vector<int> > &result, char *cache, int s, int t) {
    if (j1 >= j2 - 1) {
        if (j1 == 0) {
            result[0] = _find_i_star_list_russian(i1, j1, i2, j2, 0, seq1, seq2, cache, t, s);
        }
        if (j2 == seq2.size() - 1) {
            result[seq2.size() - 1] = _find_i_star_list_russian(i1, j1, i2, j2, seq2.size() - 1, seq1, seq2,
                                                                cache, t, s);
        }
        return;
    }

    vector<int> i_star_list = _find_i_star_list_russian(i1, j1, i2, j2, (j1 + ((j2 - j1) / 2)), seq1, seq2,
                                                        cache, t, s);

    result[(j1 + ((j2 - j1) / 2))] = i_star_list;
    _hirschberg_russians(i1, j1, i_star_list.front(), (j1 + ((j2 - j1) / 2)), seq1, seq2, scoring, result, cache, s, t);
    _hirschberg_russians(i_star_list.back(), (j1 + ((j2 - j1) / 2)), i2, j2, seq1, seq2, scoring, result, cache, s, t);
}


void print_vector(vector<int> v) {
    for (auto vi: v) {
        cout << vi << "\t";
    }
    cout << endl;
}

void print_matrix(vector<vector<int> > v) {
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[0].size(); j++) {
            cout << v[i][j] << "\t";
        }
        cout << "\n";
    }
}

pstring make_alignment(map<int, vector<int>> &result, string &seq1, string &seq2) {
    vector<pair<int, int> > flat_report;
    for (int j = 0; j < seq2.length(); j++) {
        for (auto i: result[j]) {
            flat_report.push_back(make_pair(i, j));
        }
    }

    assert(flat_report.front().first == 0);
    assert(flat_report.front().second == 0);
    vector<char> align1;
    vector<char> align2;
    for (int k = 1; k < flat_report.size(); k++) {
        pair<int, int> prev = flat_report[k - 1];
        pair<int, int> curr = flat_report[k];
        if (curr.first == prev.first) {
            assert(curr.second == prev.second + 1);
            align1.push_back('-');
            align2.push_back(seq2[curr.second]);
        } else if (curr.second == prev.second) {
            if (curr.first != prev.first + 1) {
                cerr << "Prev: " << prev.first << "\t" << prev.second << endl;
                cerr << "Curr: " << curr.first << "\t" << curr.second << endl;
            }
            assert(curr.first == prev.first + 1);
            align1.push_back(seq1[curr.first]);
            align2.push_back('-');
        } else {
            if (curr.first != prev.first + 1 || curr.second != prev.second + 1) {
                cerr << "Prev: " << prev.first << "\t" << prev.second << endl;
                cerr << "Curr: " << curr.first << "\t" << curr.second << endl;
            }
            assert(curr.first == prev.first + 1);
            assert(curr.second == prev.second + 1);
            align1.push_back(seq1[curr.first]);
            align2.push_back(seq2[curr.second]);
        }
    }
    if (flat_report.back().first != seq1.size() - 1) {
        cerr << "Problemo: " << flat_report.back().first << endl;
        cout << "Last column: ";
        print_vector(result[seq2.length() - 1]);
        cout << endl;
    }
    assert(flat_report.back().first == seq1.size() - 1);
    assert(flat_report.back().second == seq2.size() - 1);
    return make_pair(string(align1.begin(), align1.end()), string(align2.begin(), align2.end()));
}



int calculate_alignment_score(string &seq1_align, string &seq2_align, map<string, int> &scoring) {
    int score = 0;
    for (int i = 0; i < seq1_align.size(); ++i) {
        if (seq1_align[i] == '-' || seq2_align[i] == '-')
            score += scoring["gap"];
        else if (seq1_align[i] != seq2_align[i])
            score += scoring["mismatch"];
        else
            score += scoring["match"];
    }
    return score;
}

pstring hirschberg_russians(string seq1, string seq2, int s, int t , int match_score, int mismatch_score, int gap_score) {
    auto t1 = high_resolution_clock::now();
    char *cache = calculate_or_load_cache(t, s);
    auto t2 = high_resolution_clock::now();
    printf("Loading cache took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    map<string, int> scoring;
    scoring["gap"] = gap_score;
    scoring["mismatch"] = mismatch_score;
    scoring["match"] = match_score;

    if (seq1.size() < seq2.size()) // keeping sequence2 as the shorter sequence
        swap(seq1, seq2);

    map<int, vector<int> > report;
    seq1 = "0" + seq1;
    seq2 = "0" + seq2;
    _hirschberg_russians(0, 0, seq1.size() - 1, seq2.size() - 1, seq1, seq2, scoring, report, cache, s, t);

    pstring result = make_alignment(report, seq1, seq2);

    fprintf(stderr, "Our alignment score: %d\n", calculate_alignment_score(result.first, result.second, scoring));

    return result;
}


int main()
{
	string s1,s2;
	cout<<"Enter two strings\n";
    cin>>s1>>s2;
    cout<<"Enter Match, Mismatch, Gap Scores respectively\n";
    int gap_score,match_score,mismatch_score;
    cin>>match_score>>mismatch_score>>gap_score;

	Encoder e(s1 + s2);
    int s = (int) e.charset_size();
    string seq1 = e.encode(s1);
    string seq2 = e.encode(s2);

    auto t1 = high_resolution_clock::now();
    pstring alignment = hirschberg_russians(seq1, seq2, s, 2, match_score,mismatch_score,gap_score);
    auto t2 = high_resolution_clock::now();
    fprintf(stderr, "Hirschberg (Four Russians) took %lld milliseconds\n",
            duration_cast<milliseconds>(t2 - t1).count());
    cout << e.decode(alignment.first) << endl;
    cout << e.decode(alignment.second) << endl;
}