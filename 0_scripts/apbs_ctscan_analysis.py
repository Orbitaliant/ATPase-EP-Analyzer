#!/usr/bin/env python3
import os
import sys
import time
import copy
import math
import glob
import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Pool, cpu_count

# Set up logging (INFO level and above)
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

# ==============================
# FUNCTION: GET ATOMS COORDINATES AND PROTEIN DIMENSIONS
# ==============================
def prot_dim(pqr_file):
    """Get atom coordinates and protein dimensions from a PQR file (current directory)."""
    prefix = os.path.splitext(os.path.basename(pqr_file))[0]
    print(f'Getting protein dimensions for: {prefix}')
    xa, ya, za = [], [], []

    with open(pqr_file, 'r') as fr:
        lines = fr.readlines()

    dim_file = f"{prefix}_dim.dat"
    apbs_pqr_file = f"{prefix}_apbs.pqr"
    with open(dim_file, 'w') as fw, open(apbs_pqr_file, 'w') as fw2:
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                fields = line.split()
                try:
                    x = float(fields[5])
                    y = float(fields[6])
                    z = float(fields[7])
                except Exception as e:
                    print(f"Error parsing coordinates: {e}", file=sys.stderr)
                    continue
                xa.append(x)
                ya.append(y)
                za.append(z)
                fw2.write(line)

        xa = np.array(xa)
        ya = np.array(ya)
        za = np.array(za)

        if xa.size == 0 or ya.size == 0 or za.size == 0:
            raise ValueError(f"No coordinates parsed from {pqr_file}. Check file format.")

        x_len = np.max(xa) - np.min(xa)
        y_len = np.max(ya) - np.min(ya)
        z_len = np.max(za) - np.min(za)
        fo_pt = np.min(za) + z_len / 2.0

        fw.write('-PROTEIN- X len: %30.15f, Y len: %30.15f, Z len: %30.15f\n' %
                 (x_len, y_len, z_len))
        fw.write('-F0 UNIT- Z len: %30.15f, Z0: %30.15f\n' %
                 (fo_pt - np.min(za), fo_pt))

    print(f'Protein dimensions done for {prefix}.')
    return xa, ya, za, fo_pt

# ==============================
# FUNCTION: READ APBS POTENTIAL
# ==============================
def read_pot(dx_file):
    """Read APBS potential grid from a DX file in the current directory."""
    prefix = os.path.splitext(os.path.basename(dx_file))[0].replace('-pot', '')
    print(f'Reading potential grid for: {prefix}')

    with open(dx_file, 'r') as fr:
        lines = fr.readlines()

    nx = int(lines[4][35:39])
    ny = int(lines[4][39:43])
    nz = int(lines[4][43:47])
    if lines[5][7] == '-':
        ox = float(lines[5][6:20])
        if lines[5][21] == '-':
            oy = float(lines[5][20:34])
            if lines[5][35] == '-':
                oz = float(lines[5][34:48])
            else:
                oz = float(lines[5][34:47])
        else:
            oy = float(lines[5][20:33])
            if lines[5][34] == '-':
                oz = float(lines[5][33:47])
            else:
                oz = float(lines[5][33:46])
    else:
        ox = float(lines[5][6:19])
        if lines[5][20] == '-':
            oy = float(lines[5][19:33])
            if lines[5][34] == '-':
                oz = float(lines[5][33:47])
            else:
                oz = float(lines[5][33:46])
        else:
            oy = float(lines[5][19:32])
            if lines[5][33] == '-':
                oz = float(lines[5][32:46])
            else:
                oz = float(lines[5][32:45])
    dx = float(lines[6][5:18])
    dy = float(lines[7][18:31])
    dz = float(lines[8][31:44])
    xp = np.arange(ox, ox + dx*(nx - 0.5), dx)
    yp = np.arange(oy, oy + dy*(ny - 0.5), dy)
    zp = np.arange(oz, oz + dz*(nz - 0.5), dz)
    x2d = np.repeat(xp, ny).reshape((nx, ny)).T
    y2d = np.repeat(yp, nx).reshape((ny, nx))
    p = []
    total_points = nx * ny * nz
    num_full_lines = int(total_points / 3)
    for i in range(num_full_lines):
        line = lines[i + 11]
        if line[0] == '-':
            p.append(float(line[0:13]))
            if line[14] == '-':
                p.append(float(line[13:27]))
                if line[28] == '-':
                    p.append(float(line[27:41]))
                else:
                    p.append(float(line[27:40]))
            else:
                p.append(float(line[13:26]))
                if line[27] == '-':
                    p.append(float(line[26:40]))
                else:
                    p.append(float(line[26:39]))
        else:
            p.append(float(line[0:12]))
            if line[13] == '-':
                p.append(float(line[12:26]))
                if line[27] == '-':
                    p.append(float(line[26:40]))
                else:
                    p.append(float(line[26:39]))
            else:
                p.append(float(line[12:25]))
                if line[26] == '-':
                    p.append(float(line[25:39]))
                else:
                    p.append(float(line[25:38]))
    remainder = total_points % 3
    if remainder == 2:
        line = lines[num_full_lines + 11]
        if line[0] == '-':
            p.append(float(line[0:13]))
            if line[14] == '-':
                p.append(float(line[13:27]))
            else:
                p.append(float(line[13:26]))
        else:
            p.append(float(line[0:12]))
            if line[13] == '-':
                p.append(float(line[12:26]))
            else:
                p.append(float(line[12:25]))
    elif remainder == 1:
        line = lines[num_full_lines + 11]
        if line[0] == '-':
            p.append(float(line[0:13]))
        else:
            p.append(float(line[0:12]))
    p = np.array(p) * 2.57e1
    p = np.reshape(p, (nx, ny, nz))
    p = np.transpose(p)  # Ensure the first index is z

    print(f'Potential grid read for {prefix}.')
    return nx, ny, nz, xp, yp, zp, p, x2d, y2d

# ==============================
# FUNCTION: KEEP POTENTIAL WITHIN DISTANCE
# ==============================
def mask_pot(prefix, nx, ny, nz, xa, ya, za, xp, yp, zp, distmin, distmax):
    """
    Create or read a potential mask for the grid. 
    All input/output files are in the current directory.
    """
    print(f'Creating mask for: {prefix}')
    mask_path = f"{prefix}_mask.dat"
    if os.path.exists(mask_path):
        mp = []
        with open(mask_path, 'r') as fr:
            lines = fr.readlines()
        for i in range(int(nx * ny * nz / 100)):
            for j in range(100):
                mp.append(lines[i][j])
        for i in range((nz * ny * nx) % 100):
            mp.append(lines[int(nx * ny * nz / 100)][i])
        mp = np.array(mp, dtype=int).reshape(nz, ny, nx)
        mp = mp.astype(bool)
    else:
        mp = ''
        for i in range(nz):
            xa_z = xa[np.abs(za - zp[i]) <= distmax]
            ya_z = ya[np.abs(za - zp[i]) <= distmax]
            za_z = za[np.abs(za - zp[i]) <= distmax]
            if xa_z.size > 0:
                for j in range(ny):
                    condition = np.sqrt((za_z - zp[i])**2 + (ya_z - yp[j])**2) < distmax
                    xa_y = xa_z[condition]
                    ya_y = ya_z[condition]
                    za_y = za_z[condition]
                    if xa_y.size > 0:
                        for k in range(nx):
                            x_comp = xa_y - xp[k]
                            y_comp = ya_y - yp[j]
                            z_comp = za_y - zp[i]
                            comp = np.sqrt(x_comp**2 + y_comp**2 + z_comp**2)
                            comp = np.min(comp)
                            if (comp > distmin) and (comp < distmax):
                                mp = mp + '1'
                            else:
                                mp = mp + '0'
                    else:
                        mp = mp + ''.zfill(nx)
            else:
                mp = mp + ''.zfill(ny * nx)
            if ((int(100 * (i+1) / nz) != int(100 * i / nz)) and ((int(100 * (i+1) / nz) % 5) == 0)):
                print(f'{int(100 * (i+1) / nz):3d}% done for mask creation')
        print(f'Writing mask data for: {prefix}')
        with open(mask_path, 'w') as fw:
            for i in range(int(nz * ny * nx / 100)):
                fw.write(mp[i * 100:(i + 1) * 100] + '\n')
            fw.write(mp[int(nx * ny * nz / 100) * 100: int(nx * ny * nz / 100) * 100 + ((nz * ny * nx) % 100)])
        mp = np.array(list(mp), dtype=int).reshape(nz, ny, nx)
        mp = mp.astype(bool)
    print(f'Mask created for {prefix}.')
    return mp

# ==============================
# FUNCTION: PLOT POTENTIAL SURFACES ALONG Z AXIS (CT-SCAN)
# ==============================
def ct_scan(prefix, p_type, xp, yp, zp, p, x2d, y2d, mp, distmax, xa, ya, za):
    """
    Plot potential surfaces along the Z axis ("CT scan").
    Output PNG files to <p_type>_ct_scan/ in the current directory.
    """
    print(f'Plotting CT-scan surfaces for: {prefix}')
    out_dir = f"{p_type}_ct_scan"
    if os.path.exists(out_dir):
        for file in glob.glob(os.path.join(out_dir, '*.png')):
            try:
                os.remove(file)
            except Exception as e:
                print(f"Error removing file {file}: {e}")
    else:
        os.makedirs(out_dir, exist_ok=True)

    plogmax = np.log10(np.amax(np.abs(p)))
    lvl = np.concatenate((
        -1.e1 ** np.arange(np.ceil(plogmax), -5.0, -1.0),
        [0.0],
        1.e1 ** np.arange(-4.0, np.ceil(plogmax) + 0.5, 1.0)
    ), axis=None)
    lvl_contour = np.array([0])
    for i in range(len(zp)):
        fig = plt.figure()
        if (len(x2d[mp[i]]) != 0) and (len(y2d[mp[i]]) != 0):
            xmin = np.min(x2d[mp[i]]) - distmax if (np.min(x2d[mp[i]]) - distmax) > np.amin(xp) else np.amin(xp)
            xmax = np.max(x2d[mp[i]]) + distmax if (np.max(x2d[mp[i]]) + distmax) < np.amax(xp) else np.amax(xp)
            ymin = np.min(y2d[mp[i]]) - distmax if (np.min(y2d[mp[i]]) - distmax) > np.amin(yp) else np.amin(yp)
            ymax = np.max(y2d[mp[i]]) + distmax if (np.max(y2d[mp[i]]) + distmax) < np.amax(yp) else np.amax(yp)
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        cf = plt.contourf(xp, yp, p[i], levels=lvl, cmap='bwr',
                          norm=mpl.colors.SymLogNorm(linthresh=0.001))
        plt.contour(xp, yp, p[i], levels=lvl_contour, colors='black')
        fig.colorbar(cf, format='%.1e')
        dz = np.average(zp[1:] - zp[:-1])
        dclose = 0.5
        cond = (za > (zp[i] - (dclose * dz))) & (za < (zp[i] + (dclose * dz)))
        xclose = xa[cond]
        yclose = ya[cond]
        if len(xclose) != 0:
            plt.scatter(xclose, yclose, s=2, color='black')
        plt.scatter(np.average(xa), np.average(ya), s=1, color='black')
        plt.title(f'Potential (mV) CT-scan at z={zp[i]:.3f}, avg={np.average(p[i]):.3f}')
        plt.xlabel('x (Å)')
        plt.ylabel('y (Å)')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f'{i+1:05d}.png'))
        plt.close()
        # Optional progress feedback for the user
        if (10 * int(10 * (i+1) / len(zp))) != (10 * int(10 * i / len(zp))):
            print(f'{10 * int(10 * (i+1) / len(zp)):3d}% done for CT-scan of {prefix}')
    print(f'CT-scan plotting done for {prefix}.')
    return

# ==============================
# FUNCTION: CALCULATE ELECTRIC FIELD Z COMPONENT
# ==============================
def efield(prefix, zp, p):
    """
    Calculate the electric field z component from a potential grid p.
    """
    print(f'Calculating electric field z component for: {prefix}')
    ef = p.copy()
    ef[1:-1] = -(p[2:] - p[:-2]) / np.average(zp[2:] - zp[:-2])
    ef[0] = -(p[1] - p[0]) / np.average(zp[1] - zp[0])
    ef[-1] = -(p[-1] - p[-2]) / np.average(zp[-1] - zp[-2])
    print(f'Electric field calculated for {prefix}.')
    return ef

# ==============================
# FUNCTION: ANGULAR AVERAGES
# ==============================
def ang_avg(prefix, p_type, p_name, p_symb, xa, ya, x2d, y2d, zp, p, mp, fo_pt):
    """
    Calculate angular averages and plot results.
    Outputs go to <p_type>_avg/ in current directory.
    """
    print(f'Calculating angular averages for: {prefix}')
    out_dir = f"{p_type}_avg"
    if os.path.exists(out_dir):
        for file in glob.glob(os.path.join(out_dir, '*.png')):
            try:
                os.remove(file)
            except Exception as e:
                print(f"Error removing file {file}: {e}")
    else:
        os.makedirs(out_dir, exist_ok=True)
    xpc = x2d - np.average(xa)
    ypc = y2d - np.average(ya)  # Note: fixed bug (should be ya, not xa)
    zavg = []
    pavg = []
    stddev = []
    minpt = []
    maxpt = []
    nang = [20, 10]
    avg_file = os.path.join(out_dir, f"{prefix}_avg.dat")
    with open(avg_file, 'w') as fw:
        sector_file = os.path.join(out_dir, f"{prefix}_avg_sector.dat")
        with open(sector_file, 'w') as fws:
            for i in range(len(nang)):
                if i == 0:
                    print(f'Calculating fine angles for {prefix}')
                elif i == 1:
                    print(f'Calculating coarse angles for {prefix}')
                ang = np.arctan(-1.0 * xpc / ypc)
                ang = ang + (np.pi * (ypc - np.abs(ypc)) / (2.0 * ypc))
                ang = ang + (np.pi / float(nang[i]))
                ang = ang + (np.pi * (ang - np.abs(ang)) / ang)
                ang = np.array(float(nang[i]) * ang / (2.0 * np.pi), dtype=int)
                zavg_ang = [[] for _ in range(nang[i])]
                pavg_ang = [[] for _ in range(nang[i])]
                for j in range(len(zp)):
                    pts = p[j][mp[j]].copy()
                    if len(pts) != 0:
                        if i == 0:
                            zavg.append(zp[j])
                            pavg.append(np.average(pts))
                            stddev.append(np.std(pts))
                            minpt.append(np.min(pts))
                            maxpt.append(np.max(pts))
                            fw.write('%30.15e %30.15e %30.15e %30.15e %30.15e\n' %
                                     (zp[j], np.average(pts), np.std(pts), np.min(pts), np.max(pts)))
                        for k in range(nang[i]):
                            mask_ang = mp[j] * ((ang - k) == (-1.0 * (ang - k)))
                            pts_ang = p[j][mask_ang].copy()
                            if len(pts_ang) != 0:
                                zavg_ang[k].append(zp[j])
                                pavg_ang[k].append(np.average(pts_ang))
                    # Print progress every 10%
                    if (10 * int(10 * (j+1) / len(zp))) != (10 * int(10 * j / len(zp))):
                        print(f'{10 * int(10 * (j+1) / len(zp)):3d}% done for angular averages of {prefix}')
                if i == 0:  # Store sectorized data only for the finer angle resolution
                    for k in range(nang[i]):
                        if len(zavg_ang[k]) > 0:
                            fws.write(f'### Sector {k}: [{(float(k)-0.5)*360.0/float(nang[i])}°, {(float(k)+0.5)*360.0/float(nang[i])}[\n')
                            for idx in range(len(zavg_ang[k])):
                                fws.write(f'{zavg_ang[k][idx]:30.15e} {pavg_ang[k][idx]:30.15e}\n')
                            fws.write('\n')
                if i == 0:
                    zavg_np = np.array(zavg)
                    pavg_np = np.array(pavg)
                    stddev_np = np.array(stddev)
                    minpt_np = np.array(minpt)
                    maxpt_np = np.array(maxpt)
                # Instead of converting the entire list of lists into one NumPy array (which may be ragged),
                # convert each sub-list individually.
                zavg_ang = [np.array(lst) for lst in zavg_ang]
                pavg_ang = [np.array(lst) for lst in pavg_ang]
                print(f'Plotting figures for angle resolution: {nang[i]} for {prefix}')
                if i == 0:
                    for j in range(nang[i]):
                        if zavg_ang[j].size != 0:
                            plt.plot(zavg_ang[j], pavg_ang[j], label='Average', color=f'C{j % 10}')
                            plt.title(f'Average {p_name} at angles [{(float(j)-0.5)*360.0/nang[i]:.1f}°, {(float(j)+0.5)*360.0/nang[i]:.1f}°[')
                            plt.xlabel(r'$z\ (Å)$')
                            plt.ylabel(rf'${p_symb}$')
                            plt.xlim(np.min(zavg_ang[j]), np.max(zavg_ang[j]))
                            plt.legend()
                            plt.grid()
                            plt.tight_layout()
                            plt.savefig(os.path.join(out_dir, f"{prefix}_avg_ang{j:02d}.png"))
                            plt.close()
                elif i == 1:
                    for j in range(nang[i]):
                        if zavg_ang[j].size != 0:
                            plt.plot(zavg_ang[j], pavg_ang[j],
                                     label=rf'$\theta=[{(float(j)-0.5)*360.0/nang[i]:.1f}°, {(float(j)+0.5)*360.0/nang[i]:.1f}°[$',
                                     color=f'C{j % 10}')
                    plt.title(f'Average {p_name} for different angles')
                    plt.xlabel(r'$z\ (Å)$')
                    plt.ylabel(rf'${p_symb}$')
                    plt.xlim(np.min(zavg_np), np.max(zavg_np))
                    plt.legend(fontsize=8, ncol=2)
                    plt.grid()
                    plt.tight_layout()
                    plt.savefig(os.path.join(out_dir, f"{prefix}_avg_cangs.png"))
                    plt.close()
            # Plots for averages and differences
            plt.plot(zavg_np, pavg_np, label='Average')
            plt.title(f'Average {p_name}')
            plt.xlabel(r'$z\ (Å)$')
            plt.ylabel(rf'${p_symb}$')
            plt.xlim(np.min(zavg_np), np.max(zavg_np))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{prefix}_avg.png"))
            plt.close()
            plt.plot(zavg_np, pavg_np, label='Average')
            plt.plot(zavg_np, stddev_np, label='Standard Deviation')
            plt.title(f'Average {p_name} and standard deviation')
            plt.xlabel(r'$z\ (Å)$')
            plt.ylabel(rf'${p_symb}$')
            plt.xlim(np.min(zavg_np), np.max(zavg_np))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{prefix}_avg_stddev.png"))
            plt.close()
            plt.plot(zavg_np, pavg_np, label='Average')
            plt.plot(zavg_np, minpt_np, label='Minimum Point')
            plt.plot(zavg_np, maxpt_np, label='Maximum Point')
            plt.title(f'Average {p_name}, minima and maxima')
            plt.xlabel(r'$z\ (Å)$')
            plt.ylabel(rf'${p_symb}$')
            plt.xlim(np.min(zavg_np), np.max(zavg_np))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{prefix}_avg_minmax.png"))
            plt.close()
            # Differences plot
            print(f'Calculating top-mid-bot differences for {prefix}')
            midnotfound = True
            incr = 0
            while midnotfound:
                incr += 1
                if zavg_np[incr] > fo_pt:
                    avgmid = (pavg_np[incr] * (fo_pt - zavg_np[incr-1]) + pavg_np[incr-1] * (zavg_np[incr] - fo_pt)) / (zavg_np[incr] - zavg_np[incr-1])
                    midnotfound = False

            notinprot = True
            incr = 0
            delta = []
            topbot = []
            topmid = []
            midbot = []
            while notinprot:
                if zavg_np[incr] <= 0:
                    delta.append(np.abs(zavg_np[incr]))
                    topbot.append(np.abs(pavg_np[incr] - pavg_np[len(pavg_np)-1-incr]))
                    topmid.append(np.abs(pavg_np[incr] - avgmid))
                    midbot.append(np.abs(avgmid - pavg_np[len(pavg_np)-1-incr]))
                    incr += 1
                else:
                    notinprot = False
            delta = np.array(delta)
            topbot = np.array(topbot)
            topmid = np.array(topmid)
            midbot = np.array(midbot)
            plt.plot(delta, topbot, label='|Top-Bot|')
            if prefix[-2:] != 'fo':
                plt.plot(delta, topmid, label='|Top-Mid|')
                plt.plot(delta, midbot, label='|Mid-Bot|')
            plt.title(f'Average {p_name} differences')
            plt.xlabel(r'$\delta\ (Å)$')
            plt.ylabel(rf'$\Delta {p_symb}$')
            plt.xlim(0, np.max(delta))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{prefix}_avg_diff.png"))
            plt.close()
    print(f'Angular averages calculated for {prefix}.')
    return np.array(zavg), np.array(pavg)

# ==============================
# FUNCTION: Single Protein Structure Processing
# ==============================
def process_single_structure(pqr_filename, dx_filename, distmin, distmax):
    prefix = os.path.splitext(os.path.basename(pqr_filename))[0]
    print(f"Processing: {prefix}")
    xa, ya, za, fo_pt = prot_dim(pqr_filename)
    nx, ny, nz, xp, yp, zp, p, x2d, y2d = read_pot(dx_filename)
    mp = mask_pot(prefix, nx, ny, nz, xa, ya, za, xp, yp, zp, distmin, distmax)
    p_type = "pot"
    ct_scan(prefix, p_type, xp, yp, zp, p, x2d, y2d, mp, distmax, xa, ya, za)
    ef = efield(prefix, zp, p)
    p_name = "Potential"
    p_symb = "mV"
    zavg, pavg = ang_avg(prefix, p_type, p_name, p_symb, xa, ya, x2d, y2d, zp, p, mp, fo_pt)
    print(f"Overall average potential from angular averaging: {np.average(pavg):.2f}")
    return 0

# ==============================
# MAIN FUNCTION
# ==============================
def main():
    parser = argparse.ArgumentParser(description="APBS CT-scan and analysis script for a single protein (Nextflow-friendly)")
    parser.add_argument('--pqr', type=str, required=True, help="PQR filename (in current directory)")
    parser.add_argument('--dx', type=str, required=True, help="DX filename (in current directory)")
    parser.add_argument('--distmin', type=float, default=5.0, help="Minimum distance threshold (default: 5)")
    parser.add_argument('--distmax', type=float, default=15.0, help="Maximum distance threshold (default: 15)")
    args = parser.parse_args()
    process_single_structure(args.pqr, args.dx, args.distmin, args.distmax)

if __name__ == '__main__':
    main()
