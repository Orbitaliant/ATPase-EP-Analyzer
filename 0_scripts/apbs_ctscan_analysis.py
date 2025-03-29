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
def prot_dim(nam, start_time):
    prefix = os.path.basename(nam)
    logging.info('Getting protein dimensions for: %s', prefix)
    start_time2 = time.time()
    
    xa = []
    ya = []
    za = []
    
    pqr_path = os.path.join(nam, prefix + '.pqr')
    try:
        with open(pqr_path, 'r') as fr:
            lines = fr.readlines()
    except FileNotFoundError:
        logging.error("PQR file not found: %s", pqr_path)
        raise FileNotFoundError(f"PQR file not found: {pqr_path}")
    
    dim_file = os.path.join(nam, prefix + '_dim.dat')
    apbs_pqr_file = os.path.join(nam, prefix + '_apbs.pqr')
    try:
        fw = open(dim_file, 'w')
        fw2 = open(apbs_pqr_file, 'w')
    except Exception as e:
        logging.error("Error opening output files: %s", e)
        raise IOError(f"Error opening output files: {dim_file}, {apbs_pqr_file}")
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                # Use split() since your file is whitespace-delimited
                fields = line.split()
                # fields[5], fields[6], fields[7] are x, y, z
                x = float(fields[5])
                y = float(fields[6])
                z = float(fields[7])
            except Exception as e:
                logging.error("Error parsing coordinates: %s", e)
                continue
            xa.append(x)
            ya.append(y)
            za.append(z)
            fw2.write(line)
    fw2.close()
    
    xa = np.array(xa)
    ya = np.array(ya)
    za = np.array(za)
    
    if xa.size == 0 or ya.size == 0 or za.size == 0:
        logging.error("No coordinates were parsed from the PQR file %s. Check the file format.", pqr_path)
        raise ValueError(f"No coordinates parsed from {pqr_path}. Check file format.")
    
    x_len = np.max(xa) - np.min(xa)
    y_len = np.max(ya) - np.min(ya)
    z_len = np.max(za) - np.min(za)
    fo_pt = np.min(za) + z_len / 2.0
    
    fw.write('-PROTEIN- X len: %30.15f, Y len: %30.15f, Z len: %30.15f\n' %
             (x_len, y_len, z_len))
    fw.write('-F0 UNIT- Z len: %30.15f, Z0: %30.15f\n' %
             (fo_pt - np.min(za), fo_pt))
    fw.close()
    
    logging.info('Protein dimensions done for %s. Execution time: %.2f s (Total: %.2f s)',
                 prefix, time.time()-start_time2, time.time()-start_time)
    return xa, ya, za, fo_pt

# ==============================
# FUNCTION: READ APBS POTENTIAL
# ==============================
def read_pot(nam, start_time):
    prefix = os.path.basename(nam)
    logging.info('Reading potential grid for: %s', prefix)
    start_time2 = time.time()
    
    dx_path = os.path.join(nam, prefix + '-pot.dx')
    try:
        with open(dx_path, 'r') as fr:
            lines = fr.readlines()
    except FileNotFoundError:
        logging.error("DX file not found: %s", dx_path)
        raise FileNotFoundError(f"DX file not found: {dx_path}")
    
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
    
    logging.info('Potential grid read for %s. Execution time: %.2f s (Total: %.2f s)',
                 prefix, time.time()-start_time2, time.time()-start_time)
    return nx, ny, nz, xp, yp, zp, p, x2d, y2d

# ==============================
# FUNCTION: KEEP POTENTIAL WITHIN DISTANCE
# ==============================
def mask_pot(nam, nx, ny, nz, xa, ya, za, xp, yp, zp, distmin, distmax, start_time):
    prefix = os.path.basename(nam)
    logging.info('Creating mask for: %s', prefix)
    start_time2 = time.time()
    mask_path = os.path.join(nam, prefix + '_mask.dat')
    if os.path.exists(mask_path):
        mp = []
        try:
            with open(mask_path, 'r') as fr:
                lines = fr.readlines()
        except Exception as e:
            logging.error("Error reading mask file: %s", e)
            raise IOError(f"Error reading mask file: {mask_path}")
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
                logging.info('%3d%% done for mask creation', int(100 * (i+1) / nz))
        logging.info('Writing mask data for: %s', prefix)
        try:
            with open(mask_path, 'w') as fw:
                for i in range(int(nz * ny * nx / 100)):
                    fw.write(mp[i * 100:(i + 1) * 100] + '\n')
                fw.write(mp[int(nz * ny * nx / 100) * 100: int(nz * ny * nx / 100) * 100 + ((nz * ny * nx) % 100)])
        except Exception as e:
            logging.error("Error writing mask file: %s", e)
            raise IOError(f"Error writing mask file: {mask_path}")
        mp = np.array(list(mp), dtype=int).reshape(nz, ny, nx)
        mp = mp.astype(bool)
    logging.info('Mask created for %s. Execution time: %.2f s (Total: %.2f s)',
                 prefix, time.time()-start_time2, time.time()-start_time)
    return mp

# ==============================
# FUNCTION: PLOT POTENTIAL SURFACES ALONG Z AXIS (CT-SCAN)
# ==============================
def ct_scan(nam, p_type, xp, yp, zp, p, x2d, y2d, mp, distmax, xa, ya, za, start_time):
    prefix = os.path.basename(nam)
    logging.info('Plotting CT-scan surfaces for: %s', prefix)
    start_time2 = time.time()
    out_dir = os.path.join(nam, p_type + '_ct_scan')
    if os.path.exists(out_dir):
        for file in glob.glob(os.path.join(out_dir, '*.png')):
            try:
                os.remove(file)
            except Exception as e:
                logging.error("Error removing file %s: %s", file, e)
    else:
        os.makedirs(out_dir, exist_ok=True)
    
    plogmax = np.log10(np.amax(np.abs(p)))
    lvl = np.concatenate((
        -1.e1 ** np.arange(math.ceil(plogmax), -5.0, -1.0),
        [0.0],
        1.e1 ** np.arange(-4.0, math.ceil(plogmax) + 0.5, 1.0)
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
        plt.title('Potential (mV) CT-scan at z=%.3f, avg=%.3f' % (zp[i], np.average(p[i])))
        plt.xlabel('x (Å)')
        plt.ylabel('y (Å)')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, '%05d.png' % (i+1)))
        plt.close()
        if (10 * int(10 * (i+1) / len(zp))) != (10 * int(10 * i / len(zp))):
            logging.info('%3d%% done for CT-scan of %s', 10 * int(10 * (i+1) / len(zp)), prefix)
    logging.info('CT-scan plotting done for %s. Execution time: %.2f s (Total: %.2f s)',
                 prefix, time.time()-start_time2, time.time()-start_time)
    return

# ==============================
# FUNCTION: CALCULATE ELECTRIC FIELD Z COMPONENT
# ==============================
def efield(nam, zp, p, start_time):
    prefix = os.path.basename(nam)
    logging.info('Calculating electric field z component for: %s', prefix)
    start_time2 = time.time()
    ef = copy.deepcopy(p)
    ef[1:-1] = -(p[2:] - p[:-2]) / np.average(zp[2:] - zp[:-2])
    ef[0] = -(p[1] - p[0]) / np.average(zp[1] - zp[0])
    ef[-1] = -(p[-1] - p[-2]) / np.average(zp[-1] - zp[-2])
    logging.info('Electric field calculated for %s. Execution time: %.2f s (Total: %.2f s)',
                 prefix, time.time()-start_time2, time.time()-start_time)
    return ef

# ==============================
# FUNCTION: ANGULAR AVERAGES
# ==============================
def ang_avg(nam, p_type, p_name, p_symb, xa, ya, x2d, y2d, zp, p, mp, fo_pt, start_time):
    prefix = os.path.basename(nam)
    logging.info('Calculating angular averages for: %s', prefix)
    start_time2 = time.time()
    out_dir = os.path.join(nam, p_type + '_avg')
    if os.path.exists(out_dir):
        for file in glob.glob(os.path.join(out_dir, '*.png')):
            try:
                os.remove(file)
            except Exception as e:
                logging.error("Error removing file %s: %s", file, e)
    else:
        os.makedirs(out_dir, exist_ok=True)
    xpc = x2d - np.average(xa)
    ypc = y2d - np.average(xa)  # (Consider np.average(ya) if preferred)
    zavg = []
    pavg = []
    stddev = []
    minpt = []
    maxpt = []
    nang = [20, 10]
    avg_file = os.path.join(out_dir, prefix + '_avg.dat')
    try:
        fw = open(avg_file, 'w')
    except Exception as e:
        logging.error("Error opening average file: %s", e)
        raise IOError(f"Error opening average file: {avg_file}")
    sector_file = os.path.join(out_dir, prefix + '_avg_sector.dat')
    try:
        fws = open(sector_file, 'w')  # Open file for sector data
    except Exception as e:
        logging.error("Error opening sector file: %s", e)
        raise IOError(f"Error opening sector file: {sector_file}")
    for i in range(len(nang)):
        if i == 0:
            logging.info('Calculating fine angles for %s', prefix)
        elif i == 1:
            logging.info('Calculating coarse angles for %s', prefix)
        ang = np.arctan(-1.0 * xpc / ypc)
        ang = ang + (np.pi * (ypc - np.abs(ypc)) / (2.0 * ypc))
        ang = ang + (np.pi / float(nang[i]))
        ang = ang + (np.pi * (ang - np.abs(ang)) / ang)
        ang = np.array(float(nang[i]) * ang / (2.0 * np.pi), dtype=int)
        # Create lists for each angular bin (ragged lists are expected)
        zavg_ang = [[] for _ in range(nang[i])]
        pavg_ang = [[] for _ in range(nang[i])]
        for j in range(len(zp)):
            pts = copy.deepcopy(p[j][mp[j]])
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
                    pts_ang = copy.deepcopy(p[j][mask_ang])
                    if len(pts_ang) != 0:
                        zavg_ang[k].append(zp[j])
                        pavg_ang[k].append(np.average(pts_ang))
            if (10 * int(10 * (j+1) / len(zp))) != (10 * int(10 * j / len(zp))):
                logging.info('%3d%% done for angular averages of %s', 10 * int(10 * (j+1) / len(zp)), prefix)
        if i == 0:  # Store sectorized data only for the finer angle resolution
            for k in range(nang[i]):
                if len(zavg_ang[k]) > 0:
                    fws.write(f'### Sector {k}: [{(float(k)-0.5)*360.0/float(nang[i])}°, {(float(k)+0.5)*360.0/float(nang[i])}[\n')
                    for idx in range(len(zavg_ang[k])):
                        fws.write(f'{zavg_ang[k][idx]:30.15e} {pavg_ang[k][idx]:30.15e}\n')
                    fws.write('\n')
        if i == 0:
            zavg = np.array(zavg)
            pavg = np.array(pavg)
            stddev = np.array(stddev)
            minpt = np.array(minpt)
            maxpt = np.array(maxpt)
        # Instead of converting the entire list of lists into one NumPy array (which may be ragged),
        # convert each sub-list individually.
        zavg_ang = [np.array(lst) for lst in zavg_ang]
        pavg_ang = [np.array(lst) for lst in pavg_ang]
        logging.info('Plotting figures for angle resolution: %d for %s', nang[i], prefix)
        if i == 0:
            for j in range(nang[i]):
                if zavg_ang[j].size != 0:
                    plt.plot(zavg_ang[j], pavg_ang[j], label='Average', color='C%01d' % (j % 10))
                    plt.title('Average %s at angles [%.1f°, %.1f°[' %
                              (p_name, (float(j)-0.5)*360.0/float(nang[i]), (float(j)+0.5)*360.0/float(nang[i])))
                    plt.xlabel(r'$z\ (Å)$')
                    plt.ylabel(r'$%s$' % (p_symb))
                    plt.xlim(np.min(zavg_ang[j]), np.max(zavg_ang[j]))
                    plt.legend()
                    plt.grid()
                    plt.tight_layout()
                    plt.savefig(os.path.join(out_dir, prefix + '_avg_ang%02d.png' % (j)))
                    plt.close()
        elif i == 1:
            for j in range(nang[i]):
                if zavg_ang[j].size != 0:
                    plt.plot(zavg_ang[j], pavg_ang[j],
                             label=r'$\theta=[%.1f°, %.1f°[$' %
                             ((float(j)-0.5)*360.0/float(nang[i]), (float(j)+0.5)*360.0/float(nang[i])),
                             color='C%01d' % (j % 10))
            plt.title('Average %s for different angles' % (p_name))
            plt.xlabel(r'$z\ (Å)$')
            plt.ylabel(r'$%s$' % (p_symb))
            plt.xlim(np.min(zavg), np.max(zavg))
            plt.legend(fontsize=8, ncol=2)
            plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, prefix + '_avg_cangs.png'))
            plt.close()
    plt.plot(zavg, pavg, label='Average')
    plt.title('Average %s' % (p_name))
    plt.xlabel(r'$z\ (Å)$')
    plt.ylabel(r'$%s$' % (p_symb))
    plt.xlim(np.min(zavg), np.max(zavg))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, prefix + '_avg.png'))
    plt.close()
    plt.plot(zavg, pavg, label='Average')
    plt.plot(zavg, stddev, label='Standard Deviation')
    plt.title('Average %s and standard deviation' % (p_name))
    plt.xlabel(r'$z\ (Å)$')
    plt.ylabel(r'$%s$' % (p_symb))
    plt.xlim(np.min(zavg), np.max(zavg))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, prefix + '_avg_stddev.png'))
    plt.close()
    plt.plot(zavg, pavg, label='Average')
    plt.plot(zavg, minpt, label='Minimum Point')
    plt.plot(zavg, maxpt, label='Maximum Point')
    plt.title('Average %s, minima and maxima' % (p_name))
    plt.xlabel(r'$z\ (Å)$')
    plt.ylabel(r'$%s$' % (p_symb))
    plt.xlim(np.min(zavg), np.max(zavg))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, prefix + '_avg_minmax.png'))
    plt.close()
    fw.close()
    fws.close()
    logging.info('Sectorized data saved to %s', sector_file)
    
    logging.info('Calculating top-mid-bot differences for %s', prefix)
    midnotfound = True
    incr = 0
    while midnotfound:
        incr += 1
        if zavg[incr] > fo_pt:
            avgmid = (pavg[incr] * (fo_pt - zavg[incr-1]) + pavg[incr-1] * (zavg[incr] - fo_pt)) / (zavg[incr] - zavg[incr-1])
            midnotfound = False

    notinprot = True
    incr = 0
    delta = []
    topbot = []
    topmid = []
    midbot = []
    while notinprot:
        if zavg[incr] <= 0:
            delta.append(np.abs(zavg[incr]))
            topbot.append(np.abs(pavg[incr] - pavg[len(pavg)-1-incr]))
            topmid.append(np.abs(pavg[incr] - avgmid))
            midbot.append(np.abs(avgmid - pavg[len(pavg)-1-incr]))
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
    plt.title('Average %s differences' % (p_name))
    plt.xlabel(r'$\delta\ (Å)$')
    plt.ylabel(r'$\Delta %s$' % (p_symb))
    plt.xlim(0, np.max(delta))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, prefix + '_avg_diff.png'))
    plt.close()

    logging.info('Angular averages calculated for %s. Execution time: %.2f s (Total: %.2f s)',
                 prefix, time.time()-start_time2, time.time()-start_time)
    return zavg, pavg

# ==============================
# FUNCTION: Single Protein Structure Processing
# ==============================
def process_subdir(nam, start_time, distmin, distmax):
    """Processes a single subdirectory (protein structure)."""
    logging.info("Processing: %s", nam)
    
    try:
        xa, ya, za, fo_pt = prot_dim(nam, start_time)
        nx, ny, nz, xp, yp, zp, p, x2d, y2d = read_pot(nam, start_time)
        mp = mask_pot(nam, nx, ny, nz, xa, ya, za, xp, yp, zp, distmin, distmax, start_time)
        p_type = "pot"
        ct_scan(nam, p_type, xp, yp, zp, p, x2d, y2d, mp, distmax, xa, ya, za, start_time)
        ef = efield(nam, zp, p, start_time)
        p_name = "Potential"
        p_symb = "mV"
        zavg, pavg = ang_avg(nam, p_type, p_name, p_symb, xa, ya, x2d, y2d, zp, p, mp, fo_pt, start_time)

        logging.info("Processed %s", os.path.basename(nam))
        logging.info("Overall average potential from angular averaging: %.2f", np.average(pavg))
        logging.info("Execution time for %s: %.2f s", os.path.basename(nam), time.time() - start_time)

        return True  # Indicate success

    except (FileNotFoundError, IOError, ValueError) as e:
        logging.error("Error processing %s: %s", os.path.basename(nam), e)
        return False  # Indicate failure

# ==============================
# MAIN FUNCTION
# ==============================
def main():
    parser = argparse.ArgumentParser(description="APBS CT-scan and analysis script")
    # Parameterization: users can modify these values via command-line arguments
    parser.add_argument('--distmin', type=float, default=5.0, help="Minimum distance threshold (default: 5)")
    parser.add_argument('--distmax', type=float, default=15.0, help="Maximum distance threshold (default: 15)")
    parser.add_argument('--input_dir', type=str, default="4_pqr_apbs_files", help="Input directory containing PQR/APBS files")
    args = parser.parse_args()
    
    logging.info("Starting analysis with parameters: distmin=%.2f, distmax=%.2f, input_dir=%s", args.distmin, args.distmax, args.input_dir)
    
    start_time = time.time()
    
    # Get all subdirectories
    try:
        subdirs = [os.path.join(args.input_dir, d) for d in os.listdir(args.input_dir) if os.path.isdir(os.path.join(args.input_dir, d))]
    except Exception as e:
        logging.error("Error listing subdirectories in %s: %s", args.input_dir, e)
        return 1  # Return error code

    if not subdirs:
        logging.error("No subdirectories found in %s. Exiting.", args.input_dir)
        return 1  # Return error if no valid input

    num_workers = min(cpu_count(), len(subdirs))  # Limit workers to available CPU cores or subdirectory count
    logging.info("Using %d workers for parallel processing.", num_workers)

    # Prepare arguments for parallel execution
    worker_args = [(nam, start_time, args.distmin, args.distmax) for nam in subdirs]

    # Process in parallel
    with Pool(num_workers) as pool:
        results = pool.starmap(process_subdir, worker_args)  # Run tasks in parallel

    successful_runs = sum(results)  # Count successful subdirectory processing

    if successful_runs == 0:
        logging.error("All analyses failed. Exiting.")
        return 1  # Return error if no successful runs

    logging.info("All analyses completed. Total execution time: %.2f s", time.time() - start_time)
    return 0  # Return success
    
if __name__ == '__main__':
    sys.exit(main())  # Use sys.exit() to return the correct exit code
