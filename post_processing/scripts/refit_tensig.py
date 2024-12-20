import os
import subprocess

from scripts.constants import flat_results_dir


def refit_hg002_fp():
    output_path = 'tensorsignatures_data'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    base_path = os.path.join(flat_results_dir, f'STATE_HG002_vcfs/fp')
    pattern = '*.vcf$'
    process_vcf_command = PROCESS_VCF_CMD.format(base_path, pattern + ' TRUE' + f' {output_path}/hg002_fp')
    try:
        print(process_vcf_command)
        result = subprocess.check_output(process_vcf_command, shell=True, stderr=subprocess.STDOUT)
        print(result.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        # error ino
        print(f"Error: {e.output.decode('utf-8')}")
    for chunknum in ['chunk1', 'chunk2', 'chunk3']:
        tensig_prep_command = TENSIG_PREP_CMD.format(f'{output_path}/hg002_fp_{chunknum}.h5', f'{output_path}/tsdata_hg002_fp_{chunknum}.h5')
        tensig_refit_command = TENSIG_REFIT_CMD.format(f'{output_path}/tsdata_hg002_fp_{chunknum}.h5', f'{output_path}/refit_hg002_fp_{chunknum}.pkl')
        try:
            print(tensig_prep_command)
            result = subprocess.check_output(tensig_prep_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        try:
            print(tensig_refit_command)
            result = subprocess.check_output(tensig_refit_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")


def refit_subject():
    output_path = 'tensorsignatures_data'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    filternums = ['f3']
    for filternum in filternums:
        base_path = os.path.join(flat_results_dir, f'STATE_vcfs_{filternum}_region_filtered/full_genome/')
        pattern = '*.filtered3.vcf$'
        process_vcf_command = PROCESS_VCF_CMD.format(base_path, pattern + ' TRUE' +  f' {output_path}/state_{filternum}')
        tensig_prep_command = TENSIG_PREP_CMD.format(f'{output_path}/state_{filternum}.h5', f'{output_path}/tsdata_state_{filternum}.h5')
        tensig_refit_command = TENSIG_REFIT_CMD.format(f'{output_path}/tsdata_state_{filternum}.h5', f'{output_path}/refit_state_{filternum}.pkl')
        try:
            print(process_vcf_command)
            result = subprocess.check_output(process_vcf_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        for chunknum in ['chunk1', 'chunk2', 'chunk3']:
            tensig_prep_command = TENSIG_PREP_CMD.format(f'{output_path}/state_{filternum}_{chunknum}.h5', f'{output_path}/tsdata_state_{filternum}_{chunknum}.h5')
            tensig_refit_command = TENSIG_REFIT_CMD.format(f'{output_path}/tsdata_state_{filternum}_{chunknum}.h5', f'{output_path}/refit_state_{filternum}_{chunknum}.pkl')
            try:
                print(tensig_prep_command)
                result = subprocess.check_output(tensig_prep_command, shell=True, stderr=subprocess.STDOUT)
                print(result.decode('utf-8'))
            except subprocess.CalledProcessError as e:
                # error ino
                print(f"Error: {e.output.decode('utf-8')}")
            try:
                print(tensig_refit_command)
                result = subprocess.check_output(tensig_refit_command, shell=True, stderr=subprocess.STDOUT)
                print(result.decode('utf-8'))
            except subprocess.CalledProcessError as e:
                # error ino
                print(f"Error: {e.output.decode('utf-8')}")

if __name__ == "__main__":
    PROCESS_VCF_CMD = """Rscript post_processing/scripts/tensig/processVcf.R {} {}\n"""
    TENSIG_PREP_CMD = """tensorsignatures prep {} {}\n"""
    TENSIG_REFIT_CMD = """tensorsignatures refit -n {} {}\n"""
    refit_hg002_fp()
    refit_subject()



