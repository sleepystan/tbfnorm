# Copyright 2025 Stan Xiaogang Li, Dima Kozakov, Peter J. Tonge

#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.


import os
import argparse
import numpy as np
import pandas as pd
import pymol2 as pm2
from glob import glob
pd.options.mode.copy_on_write = True

# argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='directory with .pdb/.cif files, an input PDB ID, or a text file containing a list of PDB IDs. Please identify type of input!', required=True)
parser.add_argument('--input_type', '-t', help='three types of input: file_dir, pdb, pdb_list', choices=['file_dir', 'pdb', 'pdb_list'], required=True)
parser.add_argument('--out_dir', '-o', help='output directory for results')
parser.add_argument('--domain', '-d', help='Work In Progress. Functionality to truncate a protein to a specific domain in interest', default=None)
args = parser.parse_args()

def retrieve_bfactors(structures, type, out_dir, domain=None):
    with pm2.PyMOL() as s:
        if type == 'file_dir':
            files = glob(f'{args.input}/*.pdb')
            files += glob(f'{args.input}/*.cif')
            print(f'Loading {len(files)} structures...')
            
            obj_li = []
            ref = os.path.basename(files[0])[:-4]
            obj_li.append(ref)
            s.cmd.load(files[0], ref)
            for struct in files[1:]:
                obj = os.path.basename(struct)[:-4]
                if obj in obj_li:
                    raise Exception(f'Structure name {obj} is duplicated. Please check your naming!')
                s.cmd.load(struct, obj)
                s.cmd.align(obj, ref)
                obj_li.append(obj)    
        else:
            if type == 'pdb':
                s.cmd.fetch(structures)
                obj_li = [structures]
            elif type == 'pdb_list':
                pdb_li = []
                with open(structures, 'r') as fi:
                    lines = fi.readlines()
                for line in lines:
                    if not line.startswith('#'):
                        pdb_li.append(line.replace(' ', '').replace(',', '').replace('\n', ''))
                obj_li = pdb_li
                ref = pdb_li[0]
                s.cmd.fetch(ref, ref)
                for pdb in pdb_li[1:]:
                    s.cmd.fetch(pdb)
                    s.cmd.align(pdb, ref)
            else:
                raise Exception('No input recognized!')
        
        s.cmd.remove('HETATM')
        s.cmd.remove('hydrogens')
        s.cmd.remove('not (alt ""+A)')
        s.cmd.remove('not chain A')
        
        
        # additional function to remove certain residues
        # For BTK, everything except kinase domain was removed
        if False: # domain
            s.cmd.select('to_remove', 'not resi 389-659')
            
        tmp_pse_out = os.path.join(out_dir, 'tmp_structures_aln.pse')
        print(f'Saving the temporary session file to {tmp_pse_out}...')
        s.cmd.save(tmp_pse_out)
        
        # generate the first dataframe for bfactors
        ref = obj_li[0]
        bf_dict = {f'atomname': [], 'res_no': [], 'res_name': [], 'atom_id': [], f'bfactor_{ref}': []}
        s.cmd.iterate(ref, f'atomname.append(resi+resn+"_"+name)', space=bf_dict)
        s.cmd.iterate(ref, f'bfactor_{ref}.append(b)', space=bf_dict)
        s.cmd.iterate(ref, f'res_no.append(resi)', space=bf_dict)
        s.cmd.iterate(ref, f'res_name.append(resn)', space=bf_dict)
        s.cmd.iterate(ref, f'atom_id.append(name)', space=bf_dict)
        bf_df = pd.DataFrame(data=bf_dict)
        
        # generate other bfactors
        
        for obj in obj_li[1:]:
            tmp_dict = {f'atomname': [], 'res_no': [], 'res_name': [], 'atom_id': [], f'bfactor_{obj}': []}
            s.cmd.iterate(obj, f'atomname.append(resi+resn+"_"+name)', space=tmp_dict)
            s.cmd.iterate(obj, f'bfactor_{obj}.append(b)', space=tmp_dict)
            s.cmd.iterate(obj, f'res_no.append(resi)', space=tmp_dict)
            s.cmd.iterate(obj, f'res_name.append(resn)', space=tmp_dict)
            s.cmd.iterate(obj, f'atom_id.append(name)', space=tmp_dict)
            tmp_df = pd.DataFrame(data=tmp_dict)
            bf_df = bf_df.merge(tmp_df, how='outer')
            bf_df.drop_duplicates(keep='first', inplace=True)
        
        raw_out = os.path.join(out_dir, 'tmp_raw_bfactors.csv')
        print(f'Saving raw b-factors to {raw_out}...')
        bf_df.to_csv(raw_out, index=False)
        
    return bf_df, obj_li, tmp_pse_out
    
def normalize_bfactors(df, obj_li, out_dir): 
    print('Normalizing retrieved b-factors...')
    df_ = df.copy()

    means = []
    stdevs = []
    for obj in obj_li:
        tmp_array = np.array(df_[f'bfactor_{obj}'])
        avg = np.nanmean(tmp_array)
        means.append(avg)
        stdev = np.nanstd(tmp_array)
        stdevs.append(stdev)
    meanstd_dict = {'structure': obj_li, 'mean': means, 'stdev': stdevs}
    msd_df = pd.DataFrame(data=meanstd_dict)
    meanstd_out = os.path.join(out_dir, 'tmp_means_and_std.csv')
    print(f'Saving means and standard deviations to {meanstd_out}...')
    msd_df.to_csv( meanstd_out, index=False)  
    
    for i, obj in enumerate(msd_df['structure']):
        mean = msd_df.loc[i, 'mean']
        stdev = msd_df.loc[i, 'stdev']
        df_[f'normalized_{obj}'] = (df_[f'bfactor_{obj}']-mean)/stdev
        df_.drop(labels=[f'bfactor_{obj}'], axis=1, inplace=True)
    df_norm = df_
    norm_out = os.path.join(out_dir, 'normalized_bfactors.csv')
    print(f'Saving normalized b-factors to {norm_out}...')
    df_norm.to_csv(norm_out)
    
    return msd_df, df_norm

def assign_bfactors(df_norm, pse_in, obj_li, out_dir):
    pdb_out_dir = os.path.join(out_dir, 'pdbs_normalized')
    os.makedirs(pdb_out_dir, exist_ok=True)
    with pm2.PyMOL() as s:
        s.cmd.load(pse_in)
        for obj in obj_li:
            for i, resind in enumerate(df_norm['res_no']):
                resname = df_norm.loc[i, 'res_name']
                atomid = df_norm.loc[i, 'atom_id']
                if np.isnan(df_norm.loc[i, f'normalized_{obj}']):
                    # print(resind, resname, atomid)
                    continue
                s.cmd.select('atom', f'{obj} and resi {resind} and name {atomid}')
                sel_resdict = {'sel_resname': []}
                atomcount = s.cmd.iterate('atom', 'sel_resname.append(resn)', space=sel_resdict)
                # print(sel_resdict['sel_resname'][0], resind, resname, atomid)
                if sel_resdict['sel_resname'][0] != resname:
                    raise Exception(f"{sel_resdict['sel_resname'][0]} does not match {resname} in {obj}")
                if atomcount != 1:
                    raise Exception(f"More than one atom exists in {obj}: residue {resname} {resind}, atom {atomid}")
                
                b_prime = df_norm.loc[i, f"normalized_{obj}"]
                s.cmd.alter('atom', f'b={b_prime}')
            
            out_pdb = os.path.join(pdb_out_dir, f'{obj}_normalized.pdb')
            s.cmd.save(out_pdb, obj)
        # s.cmd.delete('atom')
        out_pse = os.path.join(out_dir, 'all_structures_normalized.pse')
        print(f'Saving pymol session with all structures with normalized b-factors')
        s.cmd.save(out_pse)
    return

if __name__ == '__main__':
    # Header
    print('---------------\nInitializing...\n---------------')
    statement = 'Copyright 2025 Stan Xiaogang Li, Dima Kozakov, Peter J. Tonge\n'
    statement += 'Licensed under the Apache License, Version 2.0 (the "License");\n'
    statement += 'you may not use this file except in compliance with the License.\n'
    statement += 'you may not use this file except in compliance with the License.\n'
    statement += 'You may obtain a copy of the License at\n\n'
    statement += 'http://www.apache.org/licenses/LICENSE-2.0\n\n'
    statement += 'Unless required by applicable law or agreed to in writing, software\n'
    statement += 'distributed under the License is distributed on an "AS IS" BASIS,\n'
    statement += 'WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n'
    statement += 'See the License for the specific language governing permissions and\n'
    statement += 'limitations under the License.\n'
    print(statement)

    print('_____________________________________________________\n\nIf used for a publication, please cite: [placeholder]\n\n_____________________________________________________\n\n')
    
    if not args.out_dir:
        out_dir = '.'
    elif os.path.exists(args.out_dir) and len(glob(f'{args.out_dir}/*')) != 0:
        raise Exception(f'Output exists in {args.out_dir}!')
    else:
        out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    
    raw_df, structure_list, pymol_session = retrieve_bfactors(args.input, args.input_type, out_dir)
    msd_df, df_norm = normalize_bfactors(raw_df, structure_list, out_dir)
    assign_bfactors(df_norm, pymol_session, structure_list, out_dir)
    
    # clean up fetched .cif files
    if 'pdb' in args.input_type:
        os.popen('rm ./*.cif')
    
    print('\n---------------------\nFinished calculation!\n---------------------')
    
    
        

        
        
        
            
        
        
    