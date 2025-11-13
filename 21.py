def solve(input_str: str) -> str:
    lines = input_str.strip().split('\n')
    water, time = 0, 0
    
    for i in range(1, int(lines[0]) + 1):
        T, V = map(int, lines[i].split())
        diff = T - time
        time = T
        water = water - min(water, diff) + V
    return str(water)

if __name__ == '__main__':
    data = """10
2 1
22 10
26 17
29 2
45 20
47 32
72 12
75 1
81 31
97 7
"""
    print(solve(data))