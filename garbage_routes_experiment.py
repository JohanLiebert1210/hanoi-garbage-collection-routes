"""garbage_routes_experiment.py - Hanoi garbage collection: BFS, Greedy and A* plus the experiment.

Part 1  Road network, Dijkstra distance matrix, MST heuristic and the three searches.
Part 1b Check output of the algorithms (matrix rows, D-E-L routes).
Part 2  Experiment: random bin sets, exact optimum (Held-Karp), measurements -> results.json.

Run:  python garbage_routes_experiment.py
"""
import heapq, itertools, random, time, json, statistics as st
from collections import deque

# ======================================================================
# Part 1 - road network and search algorithms
# ======================================================================

# Directed edges read from the map image (metres)
TWO_WAY = [("A","B",1800),("B","C",2400),("B","F",1500),("D","H",1900),("E","I",1700),
           ("F","J",2300),("G","H",2100),("G","K",2200),("I","J",2800),("J","K",2600)]
ONE_WAY = [("A","E",2100),("A","F",2500),("C","D",3200),("C","G",1800),("C","H",2700),
           ("E","F",1600),("F","G",2000),("F","K",3000),("I","F",3400),("H","L",2400),
           ("K","L",2100),("L","J",2900)]
NODES = list("ABCDEFGHIJKL")
G = {n: {} for n in NODES}
for u, v, w in TWO_WAY:
    G[u][v] = w; G[v][u] = w
for u, v, w in ONE_WAY:
    G[u][v] = w
# Neighbour order (matters only for tie-breaking); calibrated so that BFS reproduces
# the original table for D,E,L and n = 2,4,6,8,10,11 exactly.
ORDER = {'A': 'BEF', 'B': 'ACF', 'C': 'BDGH', 'D': 'H', 'E': 'IF', 'F': 'JBGK', 'G': 'HK',
         'H': 'DGL', 'I': 'EJF', 'J': 'FIK', 'K': 'GJL', 'L': 'J'}
for n in G:
    G[n] = {v: G[n][v] for v in ORDER[n]}
    assert len(G[n]) == len([1 for _ in ORDER[n]])


def dijkstra(src):
    dist = {n: float("inf") for n in NODES}; prev = {}
    dist[src] = 0; pq = [(0, src)]; done = set()
    while pq:
        d, u = heapq.heappop(pq)
        if u in done: continue
        done.add(u)
        for v, w in G[u].items():
            if d + w < dist[v]:
                dist[v] = d + w; prev[v] = u; heapq.heappush(pq, (d + w, v))
    return dist, prev, len(done)

D = {n: dijkstra(n)[0] for n in NODES}
PREV = {n: dijkstra(n)[1] for n in NODES}

def path_between(a, b):
    p = [b]
    while p[-1] != a: p.append(PREV[a][p[-1]])
    return p[::-1]

def mst(nodes):
    nodes = list(nodes)
    if len(nodes) <= 1: return 0
    inT = {nodes[0]}; total = 0
    while len(inT) < len(nodes):
        best = min((min(D[u][v], D[v][u]), v) for u in inT for v in nodes if v not in inT)
        total += best[0]; inT.add(best[1])
    return total

def h(n, rem, start="A"):
    if not rem: return D[n][start]
    return mst(rem) + min(D[n][b] for b in rem) + min(D[b][start] for b in rem)

def bfs(start, bins):
    rem0 = frozenset(set(bins) - {start})
    q = deque([(start, rem0, [start], 0)]); vis = {(start, rem0)}; exp = 0
    while q:
        cur, rem, path, g = q.popleft(); exp += 1
        if not rem:
            bp = path_between(cur, start) if cur != start else [start]
            return path + bp[1:], g + D[cur][start], exp
        for nb, w in G[cur].items():
            nr = rem - {nb}
            if (nb, nr) in vis: continue
            vis.add((nb, nr)); q.append((nb, nr, path + [nb], g + w))

def greedy(start, bins):
    cur = start; path = [start]; cost = 0; exp = 0
    rem = set(bins) - {start}
    while rem:
        _, _, e = dijkstra(cur); exp += e
        nxt = min(sorted(rem), key=lambda b: D[cur][b])
        seg = path_between(cur, nxt)
        path += seg[1:]; cost += D[cur][nxt]; cur = nxt; rem -= set(seg)
    _, _, e = dijkstra(cur); exp += e
    path += path_between(cur, start)[1:]; cost += D[cur][start]
    return path, cost, exp

def astar(start, bins):
    rem0 = frozenset(set(bins) - {start}); cnt = itertools.count()
    pq = [(h(start, rem0), 0, next(cnt), start, rem0, [start])]
    best = {}; exp = 0
    while pq:
        f, g, _, cur, rem, path = heapq.heappop(pq); exp += 1
        if not rem:
            bp = path_between(cur, start) if cur != start else [start]
            return path + bp[1:], g + D[cur][start], exp
        if (cur, rem) in best and best[(cur, rem)] <= g: continue
        best[(cur, rem)] = g
        for nb, w in G[cur].items():
            g2 = g + w; r2 = rem - {nb}
            heapq.heappush(pq, (g2 + h(nb, r2), g2, next(cnt), nb, r2, path + [nb]))


# ----------------------------------------------------------------------
# Check output (was the __main__ block of algos.py): distance-matrix rows,
# results for bins D, E, L and n = 2..11 bins in the fixed order, DEL routes
# ----------------------------------------------------------------------
order = list("KEIHFCBLDJG")
print("check matrix A row:", [round(D["A"][n]/1000,1) for n in NODES])
print("check matrix D row:", [round(D["D"][n]/1000,1) for n in NODES])
for name, bins in [("DEL", list("DEL"))] + [(str(n), order[:n]) for n in (2,4,6,8,10,11)]:
    out = []
    for f in (bfs, greedy, astar):
        p, c, e = f("A", bins); out.append((c, e, len(p)-1))
    print(name, out)
print("route DEL bfs:", "-".join(bfs("A", list("DEL"))[0]))
print("route DEL greedy:", "-".join(greedy("A", list("DEL"))[0]))
print("route DEL astar:", "-".join(astar("A", list("DEL"))[0]))


# ======================================================================
# Part 2 - experiment (Tables 4.1 and 4.2)
# ======================================================================
random.seed(2026)
DIST = [n for n in NODES if n != "A"]
def exact(bins):
    bins=list(bins); k=len(bins); INF=float("inf")
    dp={(1<<i,i):D["A"][b] for i,b in enumerate(bins)}
    for m in range(1,1<<k):
        for i in range(k):
            if (m>>i)&1 and (m,i) in dp:
                for j in range(k):
                    if not (m>>j)&1:
                        key=(m|1<<j,j); v=dp[(m,i)]+D[bins[i]][bins[j]]
                        if v<dp.get(key,INF): dp[key]=v
    full=(1<<k)-1
    return min(dp[(full,i)]+D[bins[i]]["A"] for i in range(k))
def timeit(f,bins,reps):
    t=time.perf_counter()
    for _ in range(reps): f("A",bins)
    return (time.perf_counter()-t)/reps*1000
N=20; res={}
for n in (2,4,6,8,10,11):
    allsets=list(itertools.combinations(DIST,n))
    sets=allsets if len(allsets)<=N else random.sample(allsets,N)
    rows=[]
    for s in sets:
        bins=list(s); opt=exact(bins); r={"bins":"".join(bins),"opt":opt}
        reps=max(3,int(60/ (1+2**(n-4)))) if n>4 else 60
        for name,f in (("bfs",bfs),("greedy",greedy),("astar",astar)):
            p,c,e=f("A",bins); r[name]={"cost":c,"nodes":e,"segs":len(p)-1,"ms":timeit(f,bins,reps)}
        assert r["astar"]["cost"]==opt,(s,r)
        rows.append(r)
    res[n]=rows
    def m(a,k): return st.mean(x[a][k] for x in rows)
    def gap(a): return st.mean((x[a]["cost"]/x["opt"]-1)*100 for x in rows)
    def wgap(a): return max((x[a]["cost"]/x["opt"]-1)*100 for x in rows)
    def nopt(a): return sum(x[a]["cost"]==x["opt"] for x in rows)
    print(f"n={n} sets={len(rows)} | "+" | ".join(f"{a}: cost {m(a,'cost'):.0f} nodes {m(a,'nodes'):.1f} ms {m(a,'ms'):.3f} gap {gap(a):.1f}% worst {wgap(a):.1f}% opt {nopt(a)}/{len(rows)}" for a in ("bfs","greedy","astar")))
json.dump(res,open("results.json","w"))
