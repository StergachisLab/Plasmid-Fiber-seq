import pandas as pd
import numpy as np
import getopt
import sys
from hmmlearn import hmm
import os
import pickle


def options():

    optlist, args = getopt.getopt(sys.argv[1:], 'i:f:e:m:',
                                  ['indir=', 'infiles=', 'context=', 'model='])

    for o, a in optlist:
        if o == '-i' or o == '--indir':
            indir = a
        elif o == '-f' or o == '--infiles':
            infiles = a.split(',')
        elif o == '-e' or o == '--context':
            context = a
        elif o == '-m' or o == '--model':
            model = a
    return indir, infiles, context, model


def encode_me(rid, read, read_info, context):
    chrom = read_info.loc[rid]['chrom']
    me = np.array(read.dropna())
    me = me[1:-1]
    # print(f'original methyl matrix = {me}')
    me = me.astype(int)
    start = read_info.loc[rid, 'start']
    end = read_info.loc[rid, 'end']
    me = me[np.where(me < (end-start))[0]]
    # print(f'Filtered methyl matrix = {me}')
    no_me = np.arange(end-start)
    no_me = np.delete(no_me, me)
    hexamers = pd.read_hdf(context, key=chrom, start=start, stop=end)
    # make sure the methylations are within the correct range
    me = me[np.where(me <= len(hexamers))[0]]
    no_me = no_me[np.where(no_me < (end-start))[0]]
    me_encode = hexamers.to_numpy().T[0]
    no_me = no_me[np.where(no_me < (me_encode.shape)[0])]
    no_me_encode = me_encode+4097
    me_encode[no_me] = 0
    me = me[np.where(me < (no_me_encode.shape)[0])]
    no_me_encode[me] = 0
    return me_encode+no_me_encode


def unpack(f_pack):
    a = np.cumsum(np.abs(f_pack))
    b = np.append(0, a)[:-1]
    c = np.zeros(max(a))+1
    d = [c[i1:i2] for i1, i2 in zip(b, a)]
    e = np.array([])
    for i in range(len(d)):
        e = np.append(e, d[i]*f_pack[i])
    return e


def apply_model(model, f, indir, outdir, context):

    read_info = pd.read_parquet(indir+f+'/'+f+'_read-info.pq')
    s_adjust = {}
    e_adjust = {}
    drops = []

    print('applying model')

    for chrom in read_info['chrom'].unique():
        print(chrom)
        fp_dic = {}
        ri = read_info.loc[read_info['chrom'] == chrom]
        reads = pd.read_parquet(indir+f+'/'+f+'_'+chrom+'.pq').T
        i = 0

        for rid, read in reads.iterrows():
            # encode the read and run it through the HMM
            read_encode = encode_me(rid, read, ri, context)
            rc = read_encode.copy()
            # check that the read encoding didn't find the wrong length (for circular)
            if len(read_encode) > 0:

                read_encode = np.append(
                    read_encode, np.append(read_encode, read_encode))

                read_encode = read_encode.astype(int).reshape(-1, 1)
                pos = model.predict(read_encode)

                # find the junctions betweens HMM states to id starts and ends of footprints
                pos_offset = np.append([0], pos)
                pos = np.append(pos, [0])
                pos_diff = pos_offset-pos
                starts = np.where(pos_diff == -1)[0]
                ends = np.where(pos_diff == 1)[0]
                ends = np.append([0], ends)

                # identify lengths of footprints and gaps based on starts/ends
                # save them as +length, -length respectively
                fps = np.sum((starts*-1, ends[1:]), axis=0).astype('int')
                gaps = -np.sum((starts, -1*ends[:-1]), axis=0).astype('int')
                combined = np.vstack((gaps, fps)).reshape((-1,), order='F')
                c = combined.copy()

                # unpacks the compressed footprints
                # splits into the middle 1/3

                # i noticed some weirdness with the last footprint not being counted correctly, so I just am grabbing the range of len(chromsome) to len(chromsome)*2
                # i hardcoded this as 16569 but you can just sub in len(pd.read_hdf(context, key=chrom)) for every 16569

                chrom_len = len(pd.read_hdf(context, key=chrom))
                if len(combined) > 0:
                    combined = unpack(combined)
                else:
                    combined = np.zeros(chrom_len)+chrom_len

                combined = combined[chrom_len:chrom_len*2]

                # this is because there were issues with reads with no footprints at all, again a slightly hacky fix. can again just sub in the chromosome length or just discard reads with NaNs after the fact
                if len(combined) < chrom_len:
                    combined = [-chrom_len]*chrom_len
                fp_dic[rid] = combined

            i += 1
            if i % 1000 == 0:
                print('completed '+str(i)+' reads')

        # export the reads as a parquet file for easy access
        # output is a parquet file where each column is a position in the circular genome, and each row is a read
        fp_df = pd.DataFrame.from_dict(fp_dic, orient='index')
        fp_df.columns = fp_df.columns.astype(str)
        print(fp_df)
        fp_df.to_parquet(outdir+f+'/'+f+'_'+chrom+'_footprints.pq')
        fp_df = ''

    # generate new read-info file to account for changes to start/end
    read_info['start'] = read_info.index.map(s_adjust)
    read_info['end'] = read_info.index.map(e_adjust)
    # remove any reads which were not long enough
    read_info = read_info.dropna(subset='start')
    read_info.to_parquet(outdir+f+'/'+f+'_read-info.pq')


indir, infiles, context, model = options()

inpq = indir+'/Infiles/parquet/'
inref = indir+'/Reference/'
context = inref+context
outdir = indir+'/Outfiles/'


with open(inref+model, 'rb') as handle:
    model = pickle.load(handle)
for f in infiles:
    print(f)
    os.system('mkdir '+outdir+f)
    apply_model(model, f, inpq, outdir, context)
