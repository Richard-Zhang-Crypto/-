import math
import time

Sbox = [0xC, 0x6, 0x9, 0x0, 0x1, 0xA, 0x2, 0xB, 0x3, 0x8, 0x5, 0xD, 0x4, 0xE, 0x7, 0xF]


def encrypt_one_round(x):
    state = [
        (x >> 28) & 0xF, (x >> 24) & 0xF, (x >> 20) & 0xF, (x >> 16) & 0xF,
        (x >> 12) & 0xF, (x >> 8) & 0xF, (x >> 4) & 0xF, (x >> 0) & 0xF
    ]

    for i in range(8):
        state[i] = Sbox[state[i]]

    t0, t1, t2, t3 = state[0], state[5], state[2], state[7]
    t4, t5, t6, t7 = state[4], state[1], state[6], state[3]

    state[0], state[1], state[2], state[3] = t0 ^ t2 ^ t3, t0, t1 ^ t2, t0 ^ t2
    state[4], state[5], state[6], state[7] = t4 ^ t6 ^ t7, t4, t5 ^ t6, t4 ^ t6

    res = 0
    res |= state[7] | (state[6] << 4) | (state[5] << 8) | (state[4] << 12)
    res |= (state[3] << 16) | (state[2] << 20) | (state[1] << 24) | (state[0] << 28)

    return res & 0xFFFFFFFF


def dot(u, y):
    z = u & y
    z ^= z >> 1
    z ^= z >> 2
    z ^= z >> 4
    z ^= z >> 8
    z ^= z >> 16
    return z & 1


def find_permutation_order():
    print("-" * 60)
    print(">>> 阶段一：计算分组密码在 2^32 状态空间下的阶 (d)")
    print("-" * 60)

    visited = bytearray(1 << 29)

    def is_visited(val):
        return (visited[val >> 3] & (1 << (val & 7))) != 0

    def set_visited(val):
        visited[val >> 3] |= (1 << (val & 7))

    global_d = 1
    cycle_count = 0
    total_space = 1 << 32
    unique_lengths = set()

    print("[过程] 开始遍历空间，追踪循环轨道 (Orbit)...")
    start_time = time.time()

    for i in range(total_space):
        if not is_visited(i):
            length = 0
            current = i

            while not is_visited(current):
                set_visited(current)
                current = encrypt_one_round(current)
                length += 1

            cycle_count += 1
            global_d = math.lcm(global_d, length)

            if length not in unique_lengths:
                unique_lengths.add(length)
                print(f"  [发现新轨道] 找到一条长度为 {length} 的循环! 当前全局阶 LCM 更新为: {global_d}")

                if length <= 100:
                    re_path = []
                    re_current = i
                    for _ in range(length):
                        re_path.append(f"0x{re_current:08X}")
                        re_current = encrypt_one_round(re_current)

                    path_str = " -> ".join(re_path)
                    print(f"  └── [轨道路径]: {path_str}")
                else:
                    print(f"  └── [轨道路径]: 长度为 {length} (> 100)，略过路径打印。")

        if i > 0 and i % 50000000 == 0:
            elapsed = time.time() - start_time
            print(f"  ... 已扫描 {i / total_space * 100:.1f}% | 耗时: {elapsed:.1f}s")

    print("\n[结果] 阶段一执行完毕！")
    print(f"  -> 总计发现独立循环个数: {cycle_count} 个")
    print(f"  -> 包含的循环长度种类有: {sorted(list(unique_lengths))}")
    print(f"  -> 最终置换的阶 (LCM) d = {global_d}")
    return global_d


def compute_all_correlations(u, v, d):
    print("\n" + "-" * 60)
    print(f">>> 阶段二：计算 1 到 {d} 轮的精确相关系数 C^(r)")
    print("-" * 60)

    counts = [0] * (d + 1)
    total_space = 1 << 32

    start_time = time.time()
    for x in range(total_space):
        y = x
        bit_u = dot(u, x)

        for r in range(1, d + 1):
            y = encrypt_one_round(y)
            bit_v = dot(v, y)
            if bit_u == bit_v:
                counts[r] += 1
            else:
                counts[r] -= 1

        if x > 0 and x % 50000000 == 0:
            elapsed = time.time() - start_time
            print(f"  ... 相关度递推计算进度: {x / total_space * 100:.1f}% | 耗时: {elapsed:.1f}s")

    cors = [0.0] * (d + 1)

    print("\n[过程] 各轮底层计数与相关度转换过程：")
    for r in range(1, d + 1):
        raw_diff = counts[r]
        match_count = (total_space + raw_diff) // 2
        mismatch_count = total_space - match_count
        cors[r] = raw_diff / float(total_space)

        print(f"  第 {r} 轮: 匹配数 {match_count} | 不匹配数 {mismatch_count}")
        print(f"         原始偏差值 = {raw_diff}")
        print(f"         相关度 C^({r}) = {raw_diff} / 2^32 = {cors[r]:.10f}")
        print("  " + "." * 40)

    return cors


def calculate_infinite_score(cors, d, lambda_weight):
    print("\n" + "-" * 60)
    print(">>> 阶段三：解析无穷轮线性壳得分 (级数求和)")
    print("-" * 60)

    print(f"[过程] 使用衰减因子 lambda = {lambda_weight}")
    print("[过程] 级数展开式: S_base = Σ (lambda^r * C^(r))")

    S_base = 0.0
    for r in range(1, d + 1):
        term = (lambda_weight ** r) * cors[r]
        S_base += term
        print(f"  -> 累加第 {r} 轮项: ({lambda_weight}^{r}) * {cors[r]:.6f} = {term:.10f}")

    print(f"\n[过程] 基础周期累加和 S_base = {S_base:.10f}")

    denominator = 1.0 - (lambda_weight ** d)
    print(f"[过程] 无穷周期缩放分母 = 1 - (lambda^{d}) = {denominator:.10f}")

    S_infinite = S_base / denominator
    print(f"[过程] 最终解析计算: S_base / 分母 = {S_infinite:.12f}")

    return S_infinite


def main():
    u = 0x000ee0f0
    v = 0x08088880
    lambda_weight = 0.5

    d = find_permutation_order()

    if d > 100:
        print(f"\n[!] 警告：置换的阶 d 过大({d})。")
        print("[!] 根据论文逻辑，此处已证明暴力算力的局限性，退出执行精确求解。")
        return

    cors = compute_all_correlations(u, v, d)
    final_score = calculate_infinite_score(cors, d, lambda_weight)

    print("\n" + "=" * 60)
    print(f"★ 结论：无穷轮线性壳解析得分 S = {final_score:.12f} ★")
    print("=" * 60)


if __name__ == "__main__":
    main()