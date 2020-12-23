import sys, os, re, math, subprocess, pickle, logging, pandas as pd
import utils
from time import sleep

# This is a wrapper for database search
# Actual database search jobs are submitted to LSF and done using databaseSearchShell.py


def generateShell(jobNumber, featureFile, paramFile, mem=1000, queue="gpu"):
    binPath = r"/hpcf/authorized_apps/proteomics_apps/jumpm/python/conda/bin/python"
    dbSearchShell = r"/hpcf/authorized_apps/proteomics_apps/jumpm/python/current/databaseSearchShell.py"
    # dbSearchShell = r"/home/jcho/dev/JUMPm/python/current/databaseSearchShell.py"
    filename = os.path.join(os.path.dirname(featureFile), "job_" + str(jobNumber) + ".sh")
    f = open(filename, "w")
    cmd = "#!/bin/bash\n"
    cmd += "#BSUB -P prot\n"
    cmd += "#BSUB -q {}\n".format(queue)
    cmd += "#BSUB -R \"rusage[mem={}]\"\n".format(mem)
    cmd += "#BSUB -e job_" + str(jobNumber) + ".e\n"
    cmd += "#BSUB -o job_" + str(jobNumber) + ".o\n"
    cmd += "{} {} {} {}\n".format(binPath, dbSearchShell, featureFile, paramFile)
    f.write(cmd)
    f.close()


def checkJobStatus(jobNumbers):
    nRemainder = len(jobNumbers)
    while nRemainder > 0:
        sleep(10)
        p = subprocess.Popen("bjobs -noheader", shell=True, stdout=subprocess.PIPE)
        (pout, _) = p.communicate()
        pout = pout.decode().strip().split("\n")
        jobRemainder = [p.split()[0] for p in pout]
        nRemainder = len(set(jobNumbers).intersection(jobRemainder))
        nFinished = len(jobNumbers) - nRemainder
        text = "\r  {} job(s) is/are finished".format(nFinished)
        sys.stdout.write(text)
        sys.stdout.flush()
    logging.info("  {} job(s) is/are finished".format(nFinished))


def submitJobs(idx, featureFile, paramFile, memSize, queue="gpu"):
    # Generation of .sh file (for each job)
    generateShell(idx, featureFile, paramFile, memSize, queue)

    # Submission of jobs to LSF
    cmd = "bsub < " + "job_" + str(idx) + ".sh"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    (pout, _) = p.communicate()
    pout = pout.decode().strip()
    pout = re.search("(?<=Job <)([0-9]+)", pout)
    jobNumber = pout.group(0)
    return jobNumber


def searchDatabase(features, paramFile, queue="gpu"):
    # features: fully-aligned features (pandas DataFrame)
    # paramFile: paramter file

    #################################
    # Preparation of job submission #
    #################################
    m = features.shape[0]
    n = 10  # Default number of entries in each job
    if int(m / n) > 200:    # When there are too many features, limit the number of jobs to 200
        n = int(m / 200) + 1
    nJobs = math.ceil(m / n)

    # Create a temporary directory for jobs (to be removed later) and change the working directory for jobs
    cwd = os.getcwd()
    tmpDir = os.path.join(cwd, ".tmp")
    if os.path.exists(tmpDir):
        os.system("rm -rf " + tmpDir)
    os.mkdir(tmpDir)
    os.system("cp " + paramFile + " " + tmpDir)    # Copy the parameter file to "tmpDir"
    os.chdir(tmpDir)    # Change the working directory to a temporary one

    ##################
    # Job submission #
    ##################
    jobNumbers = []
    mem = 1000  # Default memory reserved = 1000MB
    for i in range(nJobs):
        # Split features into "nJobs" chunks and use each chunk in each job
        start = n * i
        end = min(m, n * (i + 1))
        featureFile = "features_" + str(i) + ".pickle"
        pickle.dump(features.iloc[start:end], open(featureFile, "wb"))

        # Submission of jobs to LSF
        jobNumber = submitJobs(i, featureFile, paramFile, mem, queue)
        jobNumbers.append(jobNumber)

        text = "\r  {} job(s) is/are submitted".format(i + 1)
        sys.stdout.write(text)
        sys.stdout.flush()

    # Check the status of submitted jobs
    print()
    logging.info("  {} job(s) is/are submitted".format(nJobs))
    checkJobStatus(jobNumbers)

    ########################################################
    # Check unfinished jobs (due to not enough memory) and #
    # re-submission with the memory increase               #
    ########################################################
    print()
    print("  Checking unfinished jobs")
    logging.info("")
    logging.info("  Checking unfinished jobs")
    isFinished = False
    while not isFinished:
        jobNumbers = []
        ii = 0
        for i in range(nJobs):
            csvFile = "features_" + str(i) + ".csv"
            if not os.path.exists(csvFile): # When a job is not finished properly, there's no corresponding .csv file
                # Extracation of the required memory by parsing .o file
                f = open("job_" + str(i) + ".o")
                lines = f.read()
                mem = int(re.search("(?<=Max Memory :)\s+(\d+)", lines).group(1)) * 2    # Times 2 for safety
                f.close()

                # Re-submission of jobs
                featureFile = "features_" + str(i) + ".pickle"
                jobNumber = submitJobs(i, featureFile, paramFile, mem, queue)
                jobNumbers.append(jobNumber)

                ii += 1
                text = "\r  {} job(s) is/are submitted".format(ii + 1)
                sys.stdout.write(text)
                sys.stdout.flush()
        logging.info("  {} job(s) is/are submitted".format(ii + 1))
        # Check the status of submitted jobs
        if len(jobNumbers) > 0:
            print()
            checkJobStatus(jobNumbers)
        else:
            isFinished = True
    print()
    print("  All job(s) is/are finished")
    logging.info("")
    logging.info("  All job(s) is/are finished")

    ##########################
    # Postprocessing of jobs #
    ##########################
    res = pd.DataFrame()
    for i in range(nJobs):
        eachOutput = "features_" + str(i) + ".csv"
        try:
            df = pd.read_csv(eachOutput, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        res = res.append(df, ignore_index=True)

    ################################
    # Generation of an output file #
    ################################
    os.chdir(cwd)   # Move back to the "current working directory"
    params = utils.getParams(paramFile)
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])
    if not os.path.exists(filePath):
        os.mkdir(filePath)
    outputFile = os.path.join(filePath, "align_" + params["output_name"] + ".database_matches")
    res.to_csv(outputFile, sep="\t", index=False, na_rep="NA")

    # os.system("rm " + os.path.join(tmpDir, "features_*"))
    # os.system("rm " + os.path.join(tmpDir, "job_*"))

    return res