#define VERBOSE_ESTIMATOR 0   // 改为 1 可输出每轮搜索细节
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <string>
#include <cstdlib>
#include <chrono>
#include <cctype>

using namespace std;
using namespace std::chrono;

typedef unsigned char Cell;
typedef unsigned int uint;
typedef uint32_t mask_t;

Cell Sbox[16] = {0xC,0x6,0x9,0x0,0x1,0xA,0x2,0xB,
                 0x3,0x8,0x5,0xD,0x4,0xE,0x7,0xF};


inline int dot(uint u, uint y) {
    uint z = u & y;
    z ^= z >> 1; z ^= z >> 2; z ^= z >> 4;
    z ^= z >> 8; z ^= z >> 16;
    return z & 1;
}
inline mask_t get_nib(mask_t mask, int idx) { return (mask >> (28-4*idx)) & 0xF; }
inline mask_t set_nib(mask_t mask, int idx, mask_t val) {
    int shift = 28 - 4*idx;
    mask &= ~(0xFu << shift);
    mask |= (val & 0xF) << shift;
    return mask;
}
inline int popcnt(mask_t x) { return __builtin_popcount(x); }
inline int hdist(mask_t a, mask_t b) { return popcnt(a ^ b); }


mask_t sr(mask_t b) {
    mask_t c = 0;
    c = set_nib(c,0,get_nib(b,0)); c = set_nib(c,1,get_nib(b,5));
    c = set_nib(c,2,get_nib(b,2)); c = set_nib(c,3,get_nib(b,7));
    c = set_nib(c,4,get_nib(b,4)); c = set_nib(c,5,get_nib(b,1));
    c = set_nib(c,6,get_nib(b,6)); c = set_nib(c,7,get_nib(b,3));
    return c;
}
mask_t mc(mask_t c) {
    auto nib = [&](int i){ return get_nib(c,i); };
    mask_t d = 0;
    d = set_nib(d,0, nib(3)); d = set_nib(d,1, nib(0)^nib(1)^nib(2));
    d = set_nib(d,2, nib(1)); d = set_nib(d,3, nib(1)^nib(2)^nib(3));
    d = set_nib(d,4, nib(7)); d = set_nib(d,5, nib(4)^nib(5)^nib(6));
    d = set_nib(d,6, nib(5)); d = set_nib(d,7, nib(5)^nib(6)^nib(7));
    return d;
}
mask_t round_fwd(mask_t b) { return mc(sr(b)); }

mask_t inv_mc(mask_t a) {
    auto nib = [&](int i){ return get_nib(a,i); };
    mask_t c = 0;
    c = set_nib(c,0, nib(1)^nib(0)^nib(3)); c = set_nib(c,1, nib(2));
    c = set_nib(c,2, nib(0)^nib(2)^nib(3)); c = set_nib(c,3, nib(0));
    c = set_nib(c,4, nib(5)^nib(4)^nib(7)); c = set_nib(c,5, nib(6));
    c = set_nib(c,6, nib(4)^nib(6)^nib(7)); c = set_nib(c,7, nib(4));
    return c;
}
mask_t inv_sr(mask_t c) {
    mask_t b = 0;
    b = set_nib(b,0,get_nib(c,0)); b = set_nib(b,1,get_nib(c,5));
    b = set_nib(b,2,get_nib(c,2)); b = set_nib(b,3,get_nib(c,7));
    b = set_nib(b,4,get_nib(c,4)); b = set_nib(b,5,get_nib(c,1));
    b = set_nib(b,6,get_nib(c,6)); b = set_nib(b,7,get_nib(c,3));
    return b;
}
mask_t round_inv(mask_t a_next) { return inv_sr(inv_mc(a_next)); }


double LAT[16][16];
vector<pair<int,double>> cand_b[16], cand_a_for_b[16];

void init_LAT() {
    for (int a=0; a<16; a++)
        for (int b=0; b<16; b++) {
            int s=0;
            for (int x=0; x<16; x++) s += (dot(a,x) ^ dot(b,Sbox[x])) ? -1 : 1;
            LAT[a][b] = s / 16.0;
        }
    for (int a=0; a<16; a++) {
        for (int b=0; b<16; b++)
            if (LAT[a][b] != 0.0) cand_b[a].emplace_back(b, LAT[a][b]);
        sort(cand_b[a].begin(), cand_b[a].end(),
             [](const pair<int,double>& x, const pair<int,double>& y){ return fabs(x.second) > fabs(y.second); });
    }
    for (int b=0; b<16; b++) {
        for (int a=0; a<16; a++)
            if (LAT[a][b] != 0.0) cand_a_for_b[b].emplace_back(a, LAT[a][b]);
        sort(cand_a_for_b[b].begin(), cand_a_for_b[b].end(),
             [](const pair<int,double>& x, const pair<int,double>& y){ return fabs(x.second) > fabs(y.second); });
    }
}


inline int cand_limit(int active_cnt, int user_top_cand, int boost) {
    if (user_top_cand > 0) return user_top_cand;
    int base;
    if (active_cnt <= 4) base = 0;      
    else if (active_cnt == 5) base = 5; 
    else if (active_cnt == 6) base = 3; 
    else base = 2;                      
    if (base == 0) return 0;
    return min(base + boost, 7);
}


vector<pair<mask_t,double>> diverse_beam(unordered_map<mask_t,double>& mp, int beam) {
    vector<pair<mask_t,double>> v(mp.begin(), mp.end());
    if ((int)v.size() <= beam) {
        sort(v.begin(), v.end(), [](const pair<mask_t,double>& x, const pair<mask_t,double>& y){
            return fabs(x.second) > fabs(y.second); });
        return v;
    }
    nth_element(v.begin(), v.begin()+beam-1, v.end(),
                [](const pair<mask_t,double>& x, const pair<mask_t,double>& y){
                    return fabs(x.second) > fabs(y.second); });
    v.resize(beam);
    sort(v.begin(), v.end(), [](const pair<mask_t,double>& x, const pair<mask_t,double>& y){
        return fabs(x.second) > fabs(y.second); });

    int half = beam/2, rest = beam-half;
    vector<pair<mask_t,double>> sel;
    vector<bool> used(v.size(), false);
    for (int i=0; i<half; i++) { sel.push_back(v[i]); used[i]=true; }

    vector<mask_t> smasks; for(auto &p:sel) smasks.push_back(p.first);
    for (int k=0; k<rest; k++) {
        int best=-1, bestMinDist=-1;
        for (int i=half; i<(int)v.size(); i++) {
            if (used[i]) continue;
            int md=32;
            for (mask_t sm : smasks) md = min(md, hdist(v[i].first, sm));
            if (md > bestMinDist) { bestMinDist=md; best=i; }
        }
        if (best==-1) break;
        sel.push_back(v[best]); smasks.push_back(v[best].first); used[best]=true;
    }
    sort(sel.begin(), sel.end(), [](const pair<mask_t,double>& x, const pair<mask_t,double>& y){
        return fabs(x.second) > fabs(y.second); });
    return sel;
}


unordered_map<mask_t,double> forward(mask_t u, int rds, int max_active, int top_cand, int beam, bool silent) {
    unordered_map<mask_t,double> cur; cur[u] = 1.0;
    int boost = 0;
    for (int r=1; r<=rds; r++) {
        unordered_map<mask_t,double> nxt;
        for (auto &p : cur) {
            mask_t a = p.first; double cor = p.second;
            if (cor == 0.0) continue;
            vector<int> act; bool bad = false;
            for (int j=0; j<8; j++) {
                mask_t aj = get_nib(a,j);
                if (aj==0) continue;
                if (cand_b[aj].empty()) { bad=true; break; }
                act.push_back(j);
            }
            if (bad || (int)act.size() > max_active) continue;
            int lim = cand_limit((int)act.size(), top_cand, boost);
            vector<vector<pair<int,double>>> cands(act.size());
            for (size_t i=0; i<act.size(); i++) {
                mask_t aj = get_nib(a, act[i]);
                int cn = (lim==0) ? (int)cand_b[aj].size() : min(lim, (int)cand_b[aj].size());
                for (int k=0; k<cn; k++) cands[i].push_back(cand_b[aj][k]);
            }
            vector<int> sb(8,0);
            function<void(int,double)> dfs = [&](int idx, double prod) {
                if (idx == (int)act.size()) {
                    mask_t b = 0;
                    for (int j=0; j<8; j++) b = set_nib(b,j, sb[j]);
                    nxt[round_fwd(b)] += cor * prod;
                    return;
                }
                int pos = act[idx];
                for (auto &p_c : cands[idx]) {
                    sb[pos] = p_c.first;
                    dfs(idx+1, prod * p_c.second);
                }
            };
            dfs(0, 1.0);
        }
        auto sel = diverse_beam(nxt, beam);
        
        if ((int)sel.size() < beam * 0.3 && boost == 0) boost = 1;
        else if ((int)sel.size() > beam * 0.7 && boost == 1) boost = 0;

        if (!silent) {
            cout << "  [FWD] R" << r << ": gen=" << nxt.size() << "  retain=" << sel.size() << "  boost=" << boost << endl;
        }
        cur.clear(); for (auto &x : sel) cur[x.first] = x.second;
    }
    return cur;
}


unordered_map<mask_t,double> backward(mask_t v, int rds, int max_active, int top_cand, int beam, bool silent) {
    unordered_map<mask_t,double> cur; cur[v] = 1.0;
    int boost = 0;
    for (int r=1; r<=rds; r++) {
        unordered_map<mask_t,double> nxt;
        for (auto &p : cur) {
            mask_t anext = p.first; double cor = p.second;
            if (cor == 0.0) continue;
            mask_t b = round_inv(anext);
            vector<int> act; bool bad = false;
            for (int j=0; j<8; j++) {
                mask_t bj = get_nib(b,j);
                if (bj==0) continue;
                if (cand_a_for_b[bj].empty()) { bad=true; break; }
                act.push_back(j);
            }
            if (bad || (int)act.size() > max_active) continue;
            int lim = cand_limit((int)act.size(), top_cand, boost);
            vector<vector<pair<int,double>>> cands(act.size());
            for (size_t i=0; i<act.size(); i++) {
                mask_t bj = get_nib(b, act[i]);
                int cn = (lim==0) ? (int)cand_a_for_b[bj].size() : min(lim, (int)cand_a_for_b[bj].size());
                for (int k=0; k<cn; k++) cands[i].push_back(cand_a_for_b[bj][k]);
            }
            vector<int> sa(8,0);
            function<void(int,double)> dfs = [&](int idx, double prod) {
                if (idx == (int)act.size()) {
                    mask_t a_prev = 0;
                    for (int j=0; j<8; j++) a_prev = set_nib(a_prev,j, sa[j]);
                    nxt[a_prev] += cor * prod;
                    return;
                }
                int pos = act[idx];
                for (auto &p_c : cands[idx]) {
                    sa[pos] = p_c.first;
                    dfs(idx+1, prod * p_c.second);
                }
            };
            dfs(0, 1.0);
        }
        auto sel = diverse_beam(nxt, beam);
        if ((int)sel.size() < beam * 0.3 && boost == 0) boost = 1;
        else if ((int)sel.size() > beam * 0.7 && boost == 1) boost = 0;

        if (!silent) {
            cout << "  [BWD] R" << r << ": gen=" << nxt.size() << "  retain=" << sel.size() << "  boost=" << boost << endl;
        }
        cur.clear(); for (auto &x : sel) cur[x.first] = x.second;
    }
    return cur;
}


double estimate_correlation(int R, mask_t u, mask_t v, bool verbose) {
    int front = R/2, back = R - front;
    auto fwd = forward(u, front, 4, 0, 1000000, !verbose);   
    auto bwd = backward(v, back, 4, 0, 1000000, !verbose);
    double total = 0.0;
    for (auto &p : fwd) {
        auto it = bwd.find(p.first);
        if (it != bwd.end()) total += p.second * it->second;
    }
    return total;
}


void round_op(Cell* state) {
    for (int i=0; i<8; i++) state[i] = Sbox[state[i]];
    Cell t0=state[0], t1=state[5], t2=state[2], t3=state[7];
    Cell t4=state[4], t5=state[1], t6=state[6], t7=state[3];
    state[0]=t0^t2^t3; state[1]=t0; state[2]=t1^t2; state[3]=t0^t2;
    state[4]=t4^t6^t7; state[5]=t4; state[6]=t5^t6; state[7]=t4^t6;
}
mask_t perm(mask_t x, int R) {
    Cell state[8];
    state[0]=(x>>28)&0xF; state[1]=(x>>24)&0xF; state[2]=(x>>20)&0xF; state[3]=(x>>16)&0xF;
    state[4]=(x>>12)&0xF; state[5]=(x>>8 )&0xF; state[6]=(x>>4 )&0xF; state[7]=(x>>0 )&0xF;
    for (int r=0; r<R; r++) round_op(state);
    mask_t res=0;
    res |= (mask_t)state[7]; res |= (mask_t)state[6]<<4; res |= (mask_t)state[5]<<8; res |= (mask_t)state[4]<<12;
    res |= (mask_t)state[3]<<16; res |= (mask_t)state[2]<<20; res |= (mask_t)state[1]<<24; res |= (mask_t)state[0]<<28;
    return res;
}
double computeCor_uv(mask_t u, mask_t v, int R) {
    long long count = 0;
    const uint64_t total = 1ULL << 32;
    for (uint64_t x = 0; x < total; ++x) {
        mask_t y = perm((mask_t)x, R);
        if (dot(u, (mask_t)x) == dot(v, y)) count++;
        else count--;
    }
    return (double)count / (double)total;
}


string trim(const string &s) {
    size_t start = 0;
    while (start < s.size() && isspace((unsigned char)s[start])) start++;
    size_t end = s.size();
    while (end > start && isspace((unsigned char)s[end-1])) end--;
    return s.substr(start, end - start);
}
mask_t parse_hex(const string &raw) {
    string s = trim(raw);
    if (s.empty()) throw invalid_argument("empty hex string");
    if (s.size() >= 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))
        s = s.substr(2);
    if (s.empty()) throw invalid_argument("hex string only had prefix");
    return (mask_t)stoul(s, nullptr, 16);
}


int main() {
    init_LAT();

    string filename = "test_data.txt";
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: cannot open file '" << filename << "'" << endl;
        return 1;
    }

    ofstream fvalid("test_data_valid.txt");
    ofstream ftime("test_data_time.txt");
    ofstream fscore("test_data_score.txt");
    if (!fvalid || !ftime || !fscore) {
        cerr << "Error: cannot create output files." << endl;
        return 1;
    }

    double total_score = 0.0;
    int processed = 0, valid_cnt = 0, invalid_cnt = 0;
    string line;
    auto prog_start = high_resolution_clock::now();

    cout << "=== Batch Processor ===" << endl;
    cout << "Estimator: max_active=5, dynamic candidates + diverse beam" << endl;
    cout << "Beam width: 2,000,000" << endl;
#if VERBOSE_ESTIMATOR
    cout << "Verbose mode: ON  (showing per?round details)" << endl;
#else
    cout << "Verbose mode: OFF (set VERBOSE_ESTIMATOR to 1 for details)" << endl;
#endif
    cout << "Reading from: " << filename << endl << endl;

    while (getline(infile, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vector<string> tokens;
        {
            istringstream iss(line);
            string tok;
            while (iss >> tok) tokens.push_back(tok);
        }

        int r = 0;
        string u_str, v_str;

        if (tokens.size() == 3) {
            try { r = stoi(tokens[0]); } catch (...) { cerr << "Bad R: " << line << endl; continue; }
            u_str = tokens[1]; v_str = tokens[2];
        } else if (tokens.size() >= 4) {
            try { r = stoi(tokens[1]); } catch (...) {
                try { r = stoi(tokens[0]); u_str = tokens[1]; v_str = tokens[2]; } catch (...) { cerr << "Skipping: " << line << endl; continue; }
            }
            if (tokens.size() == 4) { u_str = tokens[2]; v_str = tokens[3]; }
        } else {
            cerr << "Skipping (expected 3+ tokens): " << line << endl;
            continue;
        }

        mask_t u, v;
        try {
            u = parse_hex(u_str);
            v = parse_hex(v_str);
        } catch (const exception &e) {
            cerr << "Hex parse error: " << line << " - " << e.what() << endl;
            continue;
        }

        cout << "[" << (processed+1) << "] R=" << r
             << "  u=0x" << hex << setw(8) << setfill('0') << u
             << "  v=0x" << hex << setw(8) << setfill('0') << v << dec << endl;

        
        auto t1 = high_resolution_clock::now();
        double V_E = estimate_correlation(r, u, v, VERBOSE_ESTIMATOR);
        auto t2 = high_resolution_clock::now();
        double ms_est = duration<double, milli>(t2 - t1).count();
        cout << "  Estimated V_E = " << scientific << setprecision(10) << V_E
             << "   (" << fixed << setprecision(3) << ms_est/1000.0 << " s)" << endl;

        
        auto t3 = high_resolution_clock::now();
        double V_T = computeCor_uv(u, v, r);
        auto t4 = high_resolution_clock::now();
        double ms_exact = duration<double, milli>(t4 - t3).count();
        cout << "  Exact    V_T = " << scientific << setprecision(10) << V_T
             << "   (" << fixed << setprecision(3) << ms_exact/1000.0 << " s)" << endl;

        processed++;

        if (fabs(V_T) > 1e-15 && fabs(V_E - V_T) < 0.25 * fabs(V_T)) {
            double score = pow(2.0, (double)r) + log2(fabs(V_E));
            total_score += score;

            fvalid << "@(" << r
                   << ",0x" << hex << setw(8) << setfill('0') << u
                   << ",0x" << hex << setw(8) << setfill('0') << v
                   << dec << "," << scientific << setprecision(10) << V_T
                   << "," << scientific << setprecision(10) << V_E << ")" << endl;

            ftime << r << " 0x" << hex << setw(8) << setfill('0') << u
                  << " 0x" << setw(8) << setfill('0') << v
                  << dec << " " << fixed << setprecision(6) << ms_est/1000.0 << " " << ms_exact/1000.0 << endl;

            fscore << scientific << setprecision(10) << score << " " << total_score << endl;

            valid_cnt++;
            cout << "  => VALID   score=" << score << "  cumulative=" << total_score << endl;
        } else {
            invalid_cnt++;
            cout << "  => INVALID (|est-exact|/|exact| = "
                 << fabs(V_E - V_T) / (fabs(V_T) + 1e-30) << ")" << endl;
        }

        auto now = high_resolution_clock::now();
        double elapsed = duration<double>(now - prog_start).count();
        cout << "  Progress: done=" << processed << "  valid=" << valid_cnt
             << "  invalid=" << invalid_cnt
             << "  elapsed=" << fixed << setprecision(1) << elapsed/3600.0 << " h" << endl << endl;
    }

    infile.close();
    fvalid.close();
    ftime.close();
    fscore.close();

    auto prog_end = high_resolution_clock::now();
    double total_time = duration<double>(prog_end - prog_start).count();

    cout << "\n========== DONE ==========" << endl;
    cout << "Total processed: " << processed << endl;
    cout << "Valid:           " << valid_cnt << endl;
    cout << "Invalid:         " << invalid_cnt << endl;
    cout << "Total score:     " << total_score << endl;
    cout << "Total time:      " << fixed << setprecision(2) << total_time/3600.0 << " hours" << endl;
    return 0;
}
