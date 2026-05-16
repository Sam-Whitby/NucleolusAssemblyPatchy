#!/usr/bin/env python3
"""
visualize_box.py

Reads a box trajectory file produced by run_box and generates a 3-panel
animated figure:
  Left:         particle positions inside the periodic square box, coloured
                by polymer type.  The box background is a uniform blue tint
                whose intensity reflects the current coupling strength γ
                (dark blue = γ≈0 denaturing; light blue = γ=1 fully coupled).
  Top-right:    energy vs simulation step.
  Bottom-right: coupling strength γ vs step (shows annealing ramp); if a stats
                file is found, accept ratio is overlaid on a twin y-axis.

Usage
-----
    python3 visualize_box.py <traj_file> [OPTIONS]

Options
    --stats FILE        Stats file (PREFIX_stats.txt).  If omitted, auto-
                        detected from the trajectory filename.
    --output FILE       Save animation to FILE (.mp4 or .gif). Default: display.
    --fps N             Frames per second (default 5).
    --skip N            Render only every N-th frame (default 1).
    --title TEXT        Figure title.

Trajectory format (written by run_box)
--------------------------------------
    <N_particles>
    step=S energy=E gamma=G L=LL nCopies=C
    <id> <poly_type> <x> <y> <copy> <ox> <oy>
    ...   (repeated N_particles times per frame)

Stats format (written by run_box)
----------------------------------
    # step  energy  acceptRatio  gamma  phase
    0  -48384.0000  0.0000  1.0000  main
    ...
"""

import argparse
import os
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np


_POLY_COLORS    = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
_BACKBONE_COLOR = "#333333"
_BOND_LW        = 1.5


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def _kv(header, key, default=None):
    m = re.search(rf'{key}=([^\s]+)', header)
    return m.group(1) if m else default


def parse_traj(path):
    """Return a list of frame dicts from a run_box trajectory file."""
    frames = []
    try:
        with open(path) as fh:
            lines = fh.readlines()
    except OSError:
        return frames

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n = int(line)
        except ValueError:
            i += 1
            continue

        i += 1
        if i >= len(lines):
            break
        hdr = lines[i].strip()
        i += 1

        step    = int(_kv(hdr, 'step',    0))
        energy  = float(_kv(hdr, 'energy', 0.0))
        gamma   = float(_kv(hdr, 'gamma',  1.0))
        L       = float(_kv(hdr, 'L',      20.0))
        nCopies = int(_kv(hdr,   'nCopies', 4))
        phase   = _kv(hdr, 'phase', 'main')

        particles = []
        for _ in range(n):
            if i >= len(lines):
                break
            parts = lines[i].strip().split()
            i += 1
            if len(parts) >= 5:
                particles.append((
                    int(parts[0]),    # id
                    int(parts[1]),    # poly_type
                    float(parts[2]),  # x
                    float(parts[3]),  # y
                    int(parts[4]),    # copy
                ))

        if len(particles) == n:
            frames.append(dict(
                step=step, energy=energy, gamma=gamma,
                L=L, nCopies=nCopies, phase=phase, particles=particles,
            ))

    return frames


def parse_stats(path):
    """Parse PREFIX_stats.txt; return list of dicts."""
    records = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        records.append(dict(
                            step=int(parts[0]),
                            energy=float(parts[1]),
                            acceptRatio=float(parts[2]),
                            gamma=float(parts[3]),
                            phase=parts[4] if len(parts) > 4 else 'main',
                        ))
                    except ValueError:
                        pass
    except OSError:
        pass
    return records


def _auto_stats(traj_path):
    """Guess PREFIX_stats.txt from PREFIX_traj.txt."""
    candidate = re.sub(r'_traj(\.txt)?$', '_stats\\1', traj_path)
    if candidate != traj_path and os.path.isfile(candidate):
        return candidate
    return None


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Animate periodic-box trajectory.")
    ap.add_argument("traj",  help="Trajectory file (output of run_box)")
    ap.add_argument("--stats",  default=None,
                    help="Stats file (default: auto-detected from traj name)")
    ap.add_argument("--output", default=None,
                    help="Save animation (.mp4 recommended; .gif is slow)")
    ap.add_argument("--fps",    type=int, default=5)
    ap.add_argument("--skip",   type=int, default=1)
    ap.add_argument("--title",  default="Periodic box assembly")
    args = ap.parse_args()

    if args.output and args.output.endswith(".gif"):
        print("Warning: GIF output can be very large and slow.\n"
              "         Consider --output file.mp4 instead (requires ffmpeg).",
              file=sys.stderr)

    print(f"Reading {args.traj} ...", flush=True)
    all_frames = parse_traj(args.traj)
    if not all_frames:
        print("No frames found.", file=sys.stderr)
        sys.exit(1)
    frames = all_frames[::args.skip]
    print(f"  {len(all_frames)} total frames, rendering {len(frames)}.", flush=True)

    stats_path = args.stats or _auto_stats(args.traj)
    stats = parse_stats(stats_path) if stats_path else []
    if stats_path and stats:
        print(f"  Stats loaded from {stats_path}.", flush=True)

    L = frames[0]['L']

    ts_steps  = [f['step']   for f in all_frames]
    ts_energy = [f['energy'] for f in all_frames]
    ts_gamma  = [f['gamma']  for f in all_frames]

    st_steps  = [s['step']        for s in stats]
    st_accept = [s['acceptRatio'] for s in stats]

    # ---------------------------------------------------------------------- #
    # Figure layout
    # ---------------------------------------------------------------------- #
    BOX_IN  = 6.0    # figure inches for the square box panel
    RIGHT_W = 2.8
    ML, MR  = 0.65, 0.30
    MT, MB  = 0.50, 0.55
    WS, HS  = 0.55, 0.45

    rp_h  = max((BOX_IN - HS) / 2.0, 0.9)
    fig_w = ML + BOX_IN + WS + RIGHT_W + MR
    fig_h = MT + BOX_IN + MB

    def fx(x): return x / fig_w
    def fy(y): return y / fig_h

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.suptitle(args.title, fontsize=10, y=1.0 - 0.1/fig_h)

    ax_anim  = fig.add_axes([fx(ML), fy(MB), fx(BOX_IN), fy(BOX_IN)])
    rx       = fx(ML + BOX_IN + WS)
    rw       = fx(RIGHT_W)
    rh       = fy(rp_h)
    ax_ener  = fig.add_axes([rx, fy(MB + rp_h + HS), rw, rh])
    ax_gamma = fig.add_axes([rx, fy(MB),              rw, rh])

    # --- Box panel static decorations ---
    pad = L * 0.02
    ax_anim.set_xlim(-pad, L + pad)
    ax_anim.set_ylim(-pad, L + pad)
    ax_anim.set_xlabel("x  (lattice units)", fontsize=9)
    ax_anim.set_ylabel("y  (lattice units)", fontsize=9)

    # Uniform blue background — intensity tracks γ.
    # Blues_r cmap: value=0 → dark blue (γ=0, denatured);
    #               value=1 → near-white (γ=1, fully coupled).
    # This matches the condensate edge colouring at full coupling.
    bg_img = ax_anim.imshow(
        np.array([[frames[0]['gamma']]]),
        extent=[0, L, 0, L],
        origin='lower', cmap='Blues_r', alpha=0.22,
        aspect='auto', zorder=0, vmin=0.0, vmax=1.0,
    )

    box_rect = mpatches.Rectangle(
        (0, 0), L, L,
        fill=False, edgecolor='steelblue', linewidth=1.5, zorder=2,
    )
    ax_anim.add_patch(box_rect)
    ax_anim.text(L + 0.4, L / 2, f"L={L:.0f}",
                 color='steelblue', fontsize=7, va='center')

    ax_anim.legend(
        handles=[mpatches.Patch(color=_POLY_COLORS[t], label=f"Polymer type {t}")
                 for t in range(4)],
        loc='upper right', fontsize=7, framealpha=0.7,
    )

    scat = ax_anim.scatter([], [], s=140, zorder=5,
                            edgecolors='k', linewidths=0.4)
    bond_lines = []
    frame_text = ax_anim.text(0.02, 0.97, "", transform=ax_anim.transAxes,
                               fontsize=8, va='top', family='monospace')

    # Equal aspect last so imshow/patches don't override it.
    ax_anim.set_aspect('equal', adjustable='box')

    # --- Energy panel ---
    ax_ener.plot(ts_steps, ts_energy, color='#555555', lw=0.8, alpha=0.4)
    ener_marker, = ax_ener.plot([], [], 'o', color='#e6194b', ms=5, zorder=5)
    ener_vline   = ax_ener.axvline(0, color='#e6194b', lw=0.8, ls='--', alpha=0.7)
    ax_ener.set_xlabel("Step", fontsize=8)
    ax_ener.set_ylabel("Energy", fontsize=8)
    ax_ener.set_title("System energy", fontsize=8)
    ax_ener.tick_params(labelsize=7)
    if ts_steps:
        ax_ener.set_xlim(min(ts_steps), max(ts_steps))

    # --- Gamma / accept-ratio panel ---
    ax_gamma.plot(ts_steps, ts_gamma, color='steelblue', lw=1.0, alpha=0.7,
                  label='γ')
    gamma_marker, = ax_gamma.plot([], [], 'o', color='steelblue', ms=5, zorder=5)
    gamma_vline   = ax_gamma.axvline(0, color='steelblue', lw=0.8, ls='--', alpha=0.7)
    ax_gamma.set_xlabel("Step", fontsize=8)
    ax_gamma.set_ylabel("γ", fontsize=8, color='steelblue')
    ax_gamma.set_title("Coupling strength γ", fontsize=8)
    ax_gamma.set_ylim(-0.05, 1.05)
    ax_gamma.tick_params(labelsize=7, colors='steelblue')
    ax_gamma.spines['left'].set_color('steelblue')
    if ts_steps:
        ax_gamma.set_xlim(min(ts_steps), max(ts_steps))

    ax2_accept = None
    accept_marker = None
    if stats:
        ax2_accept    = ax_gamma.twinx()
        ax2_accept.plot(st_steps, st_accept, color='#ff7f0e', lw=0.8, alpha=0.5,
                        label='accept ratio')
        accept_marker, = ax2_accept.plot([], [], 'o', color='#ff7f0e', ms=4, zorder=5)
        ax2_accept.set_ylabel("Accept ratio", fontsize=7, color='#ff7f0e')
        ax2_accept.tick_params(labelsize=6, colors='#ff7f0e')
        ax2_accept.spines['right'].set_color('#ff7f0e')
        ax2_accept.set_ylim(0, None)

    # ---------------------------------------------------------------------- #
    # Update function
    # ---------------------------------------------------------------------- #
    def update(frame_idx):
        nonlocal bond_lines
        fr  = frames[frame_idx]
        pts = fr['particles']
        g   = fr['gamma']
        fL  = fr['L']

        # Update uniform background intensity to match current γ.
        bg_img.set_data(np.array([[g]]))

        xs     = [p[2] for p in pts]
        ys     = [p[3] for p in pts]
        colors = [_POLY_COLORS[p[1] % 4] for p in pts]

        scat.set_offsets(np.column_stack([xs, ys]) if xs else np.empty((0, 2)))
        scat.set_color(colors)

        for ln in bond_lines:
            ln.remove()
        bond_lines = []
        by_id = {p[0]: p for p in pts}
        for p in pts:
            pid, _, px, py, copy_ = p
            if pid % 4 == 3:
                continue  # last bead of polymer; no bond to next polymer
            nxt = by_id.get(pid + 1)
            if nxt and nxt[4] == copy_:
                dx = nxt[2] - px
                dy = nxt[3] - py
                # Minimum-image for fully periodic box.
                dx -= fL * round(dx / fL)
                dy -= fL * round(dy / fL)
                if dx * dx + dy * dy < 2.5:
                    ln, = ax_anim.plot(
                        [px, px + dx], [py, py + dy],
                        color=_BACKBONE_COLOR, lw=_BOND_LW, zorder=4,
                    )
                    bond_lines.append(ln)

        frame_text.set_text(
            f"step {fr['step']}  E={fr['energy']:.0f}  γ={g:.3f}  [{fr.get('phase','main')}]"
        )

        s = fr['step']
        ener_marker.set_data([s], [fr['energy']])
        ener_vline.set_xdata([s, s])
        gamma_marker.set_data([s], [g])
        gamma_vline.set_xdata([s, s])

        artists = ([scat, frame_text, bg_img, ener_marker, ener_vline,
                    gamma_marker, gamma_vline] + bond_lines)
        if accept_marker is not None:
            # Snap accept_marker to nearest stats step.
            if st_steps:
                idx = min(range(len(st_steps)), key=lambda k: abs(st_steps[k] - s))
                accept_marker.set_data([st_steps[idx]], [st_accept[idx]])
            artists.append(accept_marker)
        return artists

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1000 // args.fps,
        blit=False,
    )

    if args.output:
        print(f"Saving to {args.output} ...", flush=True)
        writer = (animation.FFMpegWriter(fps=args.fps)
                  if args.output.endswith(".mp4")
                  else animation.PillowWriter(fps=args.fps))
        ani.save(args.output, writer=writer)
        print("Saved.")
    else:
        plt.show()


if __name__ == "__main__":
    main()
