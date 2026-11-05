# Optimal Garbage Collection Routes in Hanoi

Group 4 – Intro to AI, Autumn 2026 (Instructor: Dr. Nguyễn Trọng Nghĩa)

Members: Triệu Nguyên Bình, Lê Minh Khánh, Lê Anh Minh, Nguyễn Hà Chuyên

A garbage truck starts at the depot Ba Dinh (A), collects a set of bins in 12 Hanoi
districts (A–L) on a road network with one-way roads, and returns to A.
We compare three search algorithms:

- **Breadth-First Search (BFS)** – uninformed baseline, fewest road segments.
- **Greedy Best-First Search** – always drives to the nearest bin (Dijkstra distance).
- **A\*** – f = g + h with the admissible heuristic
  h = MST(remaining bins) + min(current → bins) + min(bins → depot).

## Run

Requires Python 3 (standard library only).

```
python garbage_routes_experiment.py
```

The script prints:

1. Check output: rows A and D of the shortest-path distance matrix, the result
   (cost, expanded nodes, road segments) of BFS, Greedy and A* for bins D, E, L and
   for n = 2, 4, 6, 8, 10, 11 bins, and the three routes for bins D, E, L.
2. Experiment: for each n, 20 random bin sets (seed 2026; all 11 sets for n = 10,
   the single set for n = 11). For every set the exact optimum is computed with a
   Held–Karp dynamic program and every A* route is checked against it. It reports
   mean cost, mean expanded nodes, mean running time, mean and worst gap to the
   optimum, and the number of sets solved optimally. Full results are saved to
   `results.json`.

Costs, node counts and gaps are deterministic; running times depend on the machine.
