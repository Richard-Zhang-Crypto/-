#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdint>
#include <string>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <numeric>
#include <sstream>
#include <climits>

using namespace std;

typedef unsigned char Cell;
typedef unsigned int uint;


Cell Sbox[16] = {0xC, 0x6, 0x9, 0x0, 0x1, 0xA, 0x2, 0xB,
                 0x3, 0x8, 0x5, 0xD, 0x4, 0xE, 0x7, 0xF};


inline int dot(uint u, uint y) {
    uint z = u & y;
    z ^= z >> 1; z ^= z >> 2; z ^= z >> 4;
    z ^= z >> 8; z ^= z >> 16;
    return z & 1;
}

inline uint get_nib(uint32_t mask, int idx) {
    return (mask >> (28 - 4 * idx)) & 0xF;
}
inline uint32_t set_nib(uint32_t mask, int idx, uint val) {
    int shift = 28 - 4 * idx;
    mask &= ~(0xFu << shift);
    mask |= (val & 0xF) << shift;
    return mask;
}


uint32_t sr_mask(uint32_t b) {
    uint nib[8];
    for (int i = 0; i < 8; ++i) nib[i] = get_nib(b, i);
    uint32_t c = 0;
    c = set_nib(c, 0, nib[0]); c = set_nib(c, 1, nib[5]);
    c = set_nib(c, 2, nib[2]); c = set_nib(c, 3, nib[7]);
    c = set_nib(c, 4, nib[4]); c = set_nib(c, 5, nib[1]);
    c = set_nib(c, 6, nib[6]); c = set_nib(c, 7, nib[3]);
    return c;
}

uint32_t mc_mask(uint32_t c) {
    uint c0 = get_nib(c,0), c1 = get_nib(c,1), c2 = get_nib(c,2), c3 = get_nib(c,3);
    uint c4 = get_nib(c,4), c5 = get_nib(c,5), c6 = get_nib(c,6), c7 = get_nib(c,7);
    uint d0 = c3, d1 = c0 ^ c1 ^ c2, d2 = c1, d3 = c1 ^ c2 ^ c3;
    uint d4 = c7, d5 = c4 ^ c5 ^ c6, d6 = c5, d7 = c5 ^ c6 ^ c7;
    uint32_t d = 0;
    d = set_nib(d, 0, d0); d = set_nib(d, 1, d1);
    d = set_nib(d, 2, d2); d = set_nib(d, 3, d3);
    d = set_nib(d, 4, d4); d = set_nib(d, 5, d5);
    d = set_nib(d, 6, d6); d = set_nib(d, 7, d7);
    return d;
}

uint32_t round_mask_prop(uint32_t b) {
    return mc_mask(sr_mask(b));
}


double LAT[16][16];
vector<pair<int, double>> cand_b[16];
vector<pair<int, double>> best_a_for_b[16];

void init_LAT() {
    for (int a = 0; a < 16; ++a)
        for (int b = 0; b < 16; ++b) {
            int sum = 0;
            for (int x = 0; x < 16; ++x)
                sum += (dot(a, x) ^ dot(b, Sbox[x])) ? -1 : 1;
            LAT[a][b] = sum / 16.0;
        }
    for (int a = 0; a < 16; ++a) {
        for (int b = 0; b < 16; ++b)
            if (LAT[a][b] != 0.0)
                cand_b[a].emplace_back(b, LAT[a][b]);
        sort(cand_b[a].begin(), cand_b[a].end(),
             [](auto &x, auto &y) { return fabs(x.second) > fabs(y.second); });
    }
    for (int b = 0; b < 16; ++b) {
        double best = 0.0;
        for (int a = 0; a < 16; ++a)
            best = max(best, fabs(LAT[a][b]));
        for (int a = 0; a < 16; ++a)
            if (fabs(LAT[a][b]) == best)
                best_a_for_b[b].emplace_back(a, LAT[a][b]);
    }
}


struct Trail {
    uint32_t u, v;
    double cor;
};

struct BeamState {
    uint32_t u;    
    uint32_t a;    
    double cor;    
};


const int TOP_K_PATHS = 1000;         


vector<Trail> gen_R1_trails() {
    vector<Trail> res;
    const int MAX_ACTIVE = 3;
    for (int act = 1; act <= MAX_ACTIVE; ++act) {
        vector<int> idx(8);
        iota(idx.begin(), idx.end(), 0);
        vector<int> chosen(act);
        function<void(int,int)> comb = [&](int start, int depth) {
            if (depth == act) {
                function<void(int,uint32_t)> assign = [&](int i, uint32_t b_val) {
                    if (i == act) {
                        uint32_t v = round_mask_prop(b_val);
                        uint32_t u = 0;
                        double cor = 1.0;
                        bool valid = true;
                        for (int j = 0; j < 8; ++j) {
                            uint bj = get_nib(b_val, j);
                            if (bj == 0) continue;
                            if (best_a_for_b[bj].empty()) { valid = false; break; }
                            auto [aj, c] = best_a_for_b[bj][0];
                            u = set_nib(u, j, aj);
                            cor *= c;
                        }
                        if (valid && u != 0 && v != 0)
                            res.push_back({u, v, cor});
                        return;
                    }
                    int pos = chosen[i];
                    for (uint bv = 1; bv < 16; ++bv)
                        assign(i+1, set_nib(b_val, pos, bv));
                };
                assign(0, 0);
            } else {
                for (int i = start; i <= 8 - act + depth; ++i) {
                    chosen[depth] = idx[i];
                    comb(i+1, depth+1);
                }
            }
        };
        comb(0, 0);
    }
    
    unordered_map<uint64_t, Trail> sum_map;
    for (auto &t : res) {
        uint64_t key = ((uint64_t)t.u << 32) | t.v;
        auto it = sum_map.find(key);
        if (it == sum_map.end())
            sum_map[key] = t;
        else
            it->second.cor += t.cor;
    }
    res.clear();
    for (auto &p : sum_map)
        res.push_back(p.second);
    return res;
}


void search_trails_beam(int R, int beam_width,
                        const vector<Trail> &init_cands,
                        vector<Trail> &results,
                        int top_cand = 0,
                        int max_active = 8)
{
    vector<BeamState> beam;
    for (auto &t : init_cands)
        beam.push_back({t.u, t.v, t.cor});
    
    sort(beam.begin(), beam.end(), [](auto &x, auto &y) {
        return fabs(x.cor) > fabs(y.cor);
    });
    if (beam_width > 0 && (int)beam.size() > beam_width)
        beam.resize(beam_width);

    for (int round = 1; round < R; ++round) {
        
        unordered_map<uint64_t, vector<double>> path_map;

        for (auto &st : beam) {
            vector<int> act_idx;
            vector<vector<pair<int,double>>> cands;
            bool invalid = false;
            for (int j = 0; j < 8; ++j) {
                uint aj = get_nib(st.a, j);
                if (aj != 0) {
                    if (cand_b[aj].empty()) { invalid = true; break; }
                    act_idx.push_back(j);
                    vector<pair<int,double>> row;
                    int cn = (top_cand == 0) ? (int)cand_b[aj].size()
                                             : min(top_cand, (int)cand_b[aj].size());
                    for (int i = 0; i < cn; ++i)
                        row.push_back(cand_b[aj][i]);
                    cands.push_back(row);
                }
            }
            if (invalid || (int)act_idx.size() > max_active) continue;

            vector<int> cur_b(8, 0);
            function<void(int,double)> dfs = [&](int idx, double cor_part) {
                if (idx == (int)act_idx.size()) {
                    uint32_t b = 0;
                    for (int j = 0; j < 8; ++j) b = set_nib(b, j, cur_b[j]);
                    uint32_t a_next = round_mask_prop(b);
                    double new_cor = st.cor * cor_part;
                    uint64_t key = ((uint64_t)st.u << 32) | a_next;
                    path_map[key].push_back(new_cor);
                    return;
                }
                int pos = act_idx[idx];
                for (auto &[bval, corval] : cands[idx]) {
                    cur_b[pos] = bval;
                    dfs(idx+1, cor_part * corval);
                }
            };
            dfs(0, 1.0);
        }

        
        unordered_map<uint64_t, BeamState> acc_map;
        for (auto &entry : path_map) {
            uint64_t key = entry.first;
            uint32_t u = (uint32_t)(key >> 32);
            uint32_t a = (uint32_t)(key & 0xFFFFFFFF);
            auto &vec = entry.second;
            
            sort(vec.begin(), vec.end(), [](double x, double y) {
                return fabs(x) > fabs(y);
            });
            
            double sum = 0.0;
            int take = min((int)vec.size(), TOP_K_PATHS);
            for (int i = 0; i < take; ++i) sum += vec[i];
            if (sum != 0.0)
                acc_map[key] = {u, a, sum};
        }

        
        beam.clear();
        for (auto &p : acc_map)
            beam.push_back(p.second);
        sort(beam.begin(), beam.end(), [](auto &x, auto &y) {
            return fabs(x.cor) > fabs(y.cor);
        });
        if (beam_width > 0 && (int)beam.size() > beam_width)
            beam.resize(beam_width);
    }

    for (auto &st : beam)
        if (st.u != 0 && st.a != 0)
            results.push_back({st.u, st.a, st.cor});
}


int main() {
    init_LAT();
    
    const int startR = 7;
    const int maxR   = 7;               
    const int beam_width = 2000000;     
    const int TOP_CAND = 4;
    const int MAX_ACTIVE = 4;

    
    vector<Trail> seeds = gen_R1_trails();
    auto sort_by_cor = [](auto &x, auto &y) { return fabs(x.cor) > fabs(y.cor); };
    sort(seeds.begin(), seeds.end(), sort_by_cor);

    int cnt_total = (int)seeds.size();
    int cnt_active[4] = {0,0,0,0};
    for (auto &t : seeds) {
        int active = 0;
        for (int i = 0; i < 8; ++i)
            if (get_nib(t.u, i) != 0) ++active;
        if (active >= 1 && active <= 3) cnt_active[active]++;
    }
    cout << "\n===== R=1 seeds statistics (max_active=3) =====" << endl;
    cout << "Total seeds : " << cnt_total << "\n\n";
    cout << "By number of active S-boxes:" << endl;
    for (int a = 1; a <= 3; ++a)
        cout << "  active = " << a << " : " << cnt_active[a] << " seeds" << endl;

    ofstream fout("linear_hulls.txt");
    fout << fixed << setprecision(10);
    fout << "# R, u(hex), v(hex), estimated_correlation\n";
    fout << "# Condition: |estimated_cor| > 1 / 2^(2^R)\n";
    fout << "# Search parameters: TOP_CAND=" << (TOP_CAND == 0 ? "all" : to_string(TOP_CAND))
         << ", MAX_ACTIVE=" << MAX_ACTIVE
         << ", BEAM_WIDTH=" << beam_width
         << ", TOP_K_PATHS=" << TOP_K_PATHS << "\n\n";

    for (int R = startR; R <= maxR; ++R) {
        vector<Trail> cands;
        if (R == 1) {
            cands = seeds;
        } else {
            search_trails_beam(R, beam_width, seeds, cands, TOP_CAND, MAX_ACTIVE);
        }

        
        unordered_map<uint64_t, Trail> uniq;
        for (auto &t : cands) {
            uint64_t key = ((uint64_t)t.u << 32) | t.v;
            auto it = uniq.find(key);
            if (it == uniq.end() || fabs(t.cor) > fabs(it->second.cor))
                uniq[key] = t;
        }

        double threshold = 1.0 / pow(2.0, pow(2.0, R));
        int passed = 0;
        for (auto &p : uniq) {
            Trail &t = p.second;
            if (t.u == 0 || t.v == 0) continue;
            if (fabs(t.cor) <= threshold) continue;

            fout << R << ", 0x" << hex << setw(8) << setfill('0') << t.u
                 << ", 0x" << setw(8) << setfill('0') << t.v
                 << dec << ", " << t.cor << endl;
            ++passed;
        }
        cout << "R=" << R << " candidates after pruning: " << passed << "/" << uniq.size() << endl;
    }

    fout.close();
    cout << "Done! Results saved to linear_hulls.txt" << endl;
    return 0;
}
