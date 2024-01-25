import pandas as pd
import numpy as np
import getopt
import argparse
import sys
from hmmlearn import hmm
import os
import pickle
import warnings
from matplotlib import pyplot as plt
warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


def options():
    # default parameters
    circular = False

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', type=str, help='input directory')
    parser.add_argument('-f', '--infiles', type=str,
                        help='input files separated by comma', default='')
    parser.add_argument('-e', '--context', type=str,
                        help='context', default='')
    parser.add_argument('-m', '--model', type=str, help='model', default='')

    args = parser.parse_args()

    indir = args.indir
    infiles = args.infiles.split(',') if args.infiles else []
    context = args.context
    model = args.model

    return indir, infiles, context, model


def encode_me(rid, read, read_info, context):

    #grab methylation info, coordinates
    chrom=read_info.loc[rid]['chrom']
    me=np.array(read.dropna())
    me=me[1:-1]
    me=me.astype(int)
    start=read_info.loc[rid,'start']
    end=read_info.loc[rid,'end']
    
    #generate array with marked methylations
    me=me[np.where(me<(end-start))[0]]
    #mask out the methylations from array for nonmethylated
    no_me=np.arange(end-start)
    no_me=np.delete(no_me, me)

    #read in hexamers for the region, encode the methylations based on the hexamers
    hexamers=pd.read_hdf(context, key=chrom, start=start, stop=end)
    me=me[np.where(me<=len(hexamers))[0]]
    no_me=no_me[np.where(no_me<(end-start))[0]]
    me_encode=hexamers.to_numpy().T[0]
    no_me = no_me[np.where(no_me < (me_encode.shape)[0])]
    no_me_encode=me_encode+4097
    #mask each other
    me_encode[no_me]=0
    me=me[np.where(me < (no_me_encode.shape)[0])]
    no_me_encode[me]=0

    #return the full array
    return me_encode+no_me_encode


def unpack(f_pack):
    if len(f_pack) == 0:
        return np.array([])

    a = np.cumsum(np.abs(f_pack))
    b = np.append(0, a)[:-1]
    c = np.zeros(max(a))+1
    d = [c[i1:i2] for i1, i2 in zip(b, a)]
    e = np.array([])
    for i in range(len(d)):
        e = np.append(e, d[i]*f_pack[i])
    return e

# if not a.any():
#         print(f'Found an empty array -> {a} - Continuing...')
#         return


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
        # read in hexamers for the chromosome
        me_encode = pd.read_hdf(context, key=chrom).to_numpy().T[0]

        for rid, read in reads.iterrows():
            # encode the read and run it through the HMM
            read_encode = encode_me(rid, read, ri, me_encode)
            f'encode_me output = {print(np.unique(read_encode))}'
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

                # unpacks the compressed footprints
                # splits into the middle 1/3
                combined = unpack(combined)
                combined = combined[len(combined)//3:(len(combined)*2//3)]
                new_s = 0
                new_e = len(me_encode)

                s_adjust[rid] = new_s
                e_adjust[rid] = new_e
                fp_dic[rid] = combined
            else:
                s_adjust[rid] = np.nan
                e_adjust[rid] = np.nan

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
