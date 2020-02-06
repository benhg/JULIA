def run_single_align(filename):
    output = 0

    uniquely = os.popen(
        "awk 'FNR == 9{print $6}' " +
        strsvalue +
        "Log.final.out").read()
    if uniquely == 0:
        uniquely = '1'
    uniquely = int(uniquely)

    multi = os.popen(
        "awk 'FNR == 24{print $9}' " +
        strsvalue +
        "Log.final.out").read()
    if multi == '0':
        multi = '1'
    multi = int(multi)

    totalreads = os.popen(
        "awk 'FNR == 6{print $6}' " +
        strsvalue +
        "Log.final.out").read()

    if uniquely + multi >= threshold:
        os.system("echo " + filename + f">> {directory}/final.txt")
        output = 1
        uniquelyC = 'NO'
        multiC = 'NO'
        totalreadsC = 'NO'
        percRatio = 'NO'

    if svalue < 12:
        counter = 1
    else:
        counter = 15

    while output == 0 and counter > 0:
        os.system("echo" + filename + f" >> {directory}/maybehopper.txt")
        if counter == svalue:
            counter += 1
        if counter == 20 or counter == 7:
            counter = 0
        elif counter == 17:
            counter = 19

        strcounter = str(counter)

        star_align(directory, mvalue, strcounter, genomeDir)

        threshold = 1000

        uniquelyC = os.popen(
            "awk 'FNR == 9{print $6}' " +
            strcounter +
            "Log.final.out").read()
        if uniquelyC == '0':
            uniquelyC = '1'
        uniquelyC = int(uniquelyC)

        multiC = os.popen(
            "awk 'FNR == 24{print $9}' " +
            strcounter +
            "Log.final.out").read()
        if multiC == '0':
            multiC = '1'
        multiC = int(multiC)

        totalreadsC = os.popen(
            "awk 'FNR == 6{print $6}' " +
            strcounter +
            "Log.final.out").read()

        percRatio = (multiC + uniquelyC) / (multi + uniquely)

        if (multiC + uniquelyC) > threshold and percRatio > 100:
            output = 1
            os.system("echo " + filename + f" >> {directory}/hopper.txt")
        else:
            if counter == 17:
                counter = 19
            elif counter == 20 or counter == 7:
                counter = 0
            else:
                counter = counter + 1

    if output == 1:  # it needs to be added to a file
        os.system(
            "echo " +
            str(filename) +
            "  " +
            str(uniquely) +
            "  " +
            str(multi) +
            "  " +
            str(totalreads) +
            "  " +
            str(uniquelyC) +
            "  " +
            str(multiC) +
            "  " +
            str(totalreadsC) +
            "  " +
            str(percRatio) +
            f" >> {directory}/index_hopping_output.txt")

    if output == 0:
        os.system("echo " + filename + f" >> {directory}/never.txt")