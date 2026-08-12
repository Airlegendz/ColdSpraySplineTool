"""Small matplotlib helper shared by nozzle_geometry.py and sweep.py."""

from typing import List, Tuple

Point = Tuple[float, float]


def plot_profile(points: List[Point], out_path: str, title: str = "") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = [p[0] for p in points]
    rs = [p[1] for p in points]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(xs, rs, "b-", linewidth=1.5)
    ax.plot(xs, [-r for r in rs], "b-", linewidth=1.5)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("r (m)")
    ax.set_aspect("equal", adjustable="datalim")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
