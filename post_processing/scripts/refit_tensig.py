import os
import subprocess
import sys

from constants import flat_results_dir


def refit_hg002_fp():
    PROCESS_VCF_CMD = """Rscript tensig/processVcf.R {} {}\n"""
    TENSIG_PREP_CMD = """tensorsignatures prep {} {}\n"""
    TENSIG_REFIT_CMD = """tensorsignatures refit -n {} {}\n"""
    output_path = 'tensorsignatures_data'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    base_path = os.path.join(flat_results_dir, f'STATE_HG002_vcfs/fp')
    pattern = '*.vcf$'
    process_vcf_command = PROCESS_VCF_CMD.format(base_path, pattern + ' FALSE' + f' {output_path}/hg002_fp')
    tensig_prep_command = TENSIG_PREP_CMD.format(f'{output_path}/hg002_fp.h5', f'{output_path}/tsdata_hg002_fp.h5')
    tensig_refit_command = TENSIG_REFIT_CMD.format(f'{output_path}/tsdata_hg002_fp.h5', f'{output_path}/refit_hg002_fp.pkl')
    try:
        print(process_vcf_command)
        result = subprocess.check_output(process_vcf_command, shell=True, stderr=subprocess.STDOUT)
        print(result.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        # error ino
        print(f"Error: {e.output.decode('utf-8')}")
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
    latest_run = str(sys.argv[1])
    PROCESS_VCF_CMD = """Rscript tensig/processVcf.R {} {}\n"""
    TENSIG_PREP_CMD = """tensorsignatures prep {} {}\n"""
    TENSIG_REFIT_CMD = """tensorsignatures refit -n {} {}\n"""
    output_path = 'tensorsignatures_data'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    filternums = ['f1', 'f3']
    for filternum in filternums:
        base_path = os.path.join(flat_results_dir, f'STATE_vcfs_{filternum}_region_filtered')
        pattern = '*.filtered[1,3].vcf$'
        process_vcf_command = PROCESS_VCF_CMD.format(base_path, pattern + ' FALSE' +  f' {output_path}/state{latest_run}_{filternum}')
        tensig_prep_command = TENSIG_PREP_CMD.format(f'{output_path}/state{latest_run}_{filternum}.h5', f'{output_path}/tsdata_state{latest_run}_{filternum}.h5')
        tensig_refit_command = TENSIG_REFIT_CMD.format(f'{output_path}/tsdata_state{latest_run}_{filternum}.h5', f'{output_path}/refit_state{latest_run}_{filternum}.pkl')
        try:
            print(process_vcf_command)
            result = subprocess.check_output(process_vcf_command, shell=True, stderr=subprocess.STDOUT)
            print(result.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            # error ino
            print(f"Error: {e.output.decode('utf-8')}")
        for chunknum in ['chunk1', 'chunk2', 'chunk3']:
            tensig_prep_command = TENSIG_PREP_CMD.format(f'{output_path}/state{latest_run}_{filternum}_{chunknum}.h5', f'{output_path}/tsdata_state{latest_run}_{filternum}_{chunknum}.h5')
            tensig_refit_command = TENSIG_REFIT_CMD.format(f'{output_path}/tsdata_state{latest_run}_{filternum}_{chunknum}.h5', f'{output_path}/refit_state{latest_run}_{filternum}_{chunknum}.pkl')
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
    refit_hg002_fp()
    refit_subject()




# input_vcf_paths1 = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', '*', '*.ann.filtered1.vcf.gz'))
# input_vcf_paths1.sort()
# input_list1 = ' '.join(input_vcf_paths1)

# input_vcf_paths3 = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', '*', '*.ann.filtered3.vcf.gz'))
# input_vcf_paths3.sort()
# input_list3 = ' '.join(input_vcf_paths3)


# input_fp_vcf_paths = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', 'HG002_*', 'fp', '0000.vcf'))
# input_fp_vcf_paths.sort()
# input_fp_list = ' '.join(input_fp_vcf_paths)


# input_tn_vcf_paths = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', 'HG002_*', 'fp', '0001.vcf'))
# input_tn_vcf_paths.sort()



# input_tp_vcf_paths = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', 'HG002_*', 'fp', '0002.vcf'))
# input_tp_vcf_paths.sort()
# input_tpfptn_list = ' '.join(input_tp_vcf_paths + input_fp_vcf_paths + input_tn_vcf_paths)

# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_tpfptn_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_tpfptn_list, 'hg002_tpfptn.h5'))



# input_fp_pileup_paths = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', 'HG002_*', 'fp_pileup', '0000.vcf'))
# input_fp_pileup_paths.sort()
# input_fp_pileup_list = ' '.join(input_fp_pileup_paths)





# input_vcf_fail_paths1 = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', '*', '*.ann.failed_filter1.vcf.gz'))
# input_vcf_fail_paths1.sort()

# input_vcf_fail_paths3 = glob.glob(os.path.join(base_path, '*', 'results/variant_calling/clairs', '*', '*.ann.failed_filter3.vcf.gz'))
# input_vcf_fail_paths3.sort()
# input_vcf_fail_list1 = ' '.join(input_vcf_fail_paths1)
# input_vcf_fail_list3 = ' '.join( input_vcf_fail_paths3)

# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_patient_vcf_f1.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_list1, 'state002_{latest_run}_f1.h5'))

# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_patient_vcf_{filternum}.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_list3, 'state002_{latest_run}_{filternum}.h5'))


# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_patient_failed_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_vcf_fail_list1, 'state002_032_failed.h5'))



# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_fp_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_fp_list, 'hg002_fp.h5'))


# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_fp_pileup_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_fp_pileup_list, 'hg002_fp_pileup_list.h5'))


# input_patient_fp_list = ' '.join(input_vcf_paths + input_fp_vcf_paths)
# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_patients_fp_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_patient_fp_list, 'state002_029_fp.h5'))



# input_patient_fp_fp_pileup_list = ' '.join(input_vcf_paths + input_fp_vcf_paths + input_fp_pileup_paths)
# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_patients_fp_fp_pileupe_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_patient_fp_fp_pileup_list, 'state002_029_fp_fp_pileup.h5'))




# input_patient_patient_fail_list = ' '.join(input_vcf_paths + input_vcf_fail_paths)
# with open(os.path.join('/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/process_patients_filtered_failed_vcf.sh'), 'w') as bash_file:
#     bash_file.write(PROCESS_VCF_CMD.format(input_patient_patient_fail_list, 'state002_029_filtered_failed.h5'))

