### SCRIPT TO EXTRACT THE INTRONIC PERCENT FROM QUALIMAP REPORT ###

with open("intronic_info.txt") as f:
    with open("intronic_percent.txt", "a") as f2:
        for line in f:
            percent = line.split(" ")[-1]
            percent = percent[1:]
            percent = percent[0:-3]
            print(percent)
            f2.write(percent)
            f2.write("\n")
