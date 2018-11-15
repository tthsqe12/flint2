"""
Make a "broken" horizontal bar plot, i.e., one with gaps
"""
import matplotlib.pyplot as plt

data = [[0,[0,0,1,0,2,0,3,8,4,8,5,27,6,27,7,28,8,28,9,62,10,62,11,63,12,63,13,121,14,121,15,121,16,121,17,134,18,134,19,134,20,134,21,281,22,281,23,282,24,282,25,443,26,443,27,443,28,443,29,630,30,630,31,630,32,630,33,691,34,691,35,691,36,691,37,691,38,691,39,691,40,691,41,832,42,832,43,832,44,832,45,1054,46,1054,47,1054,48,1054,49,1311,50,1311,51,1311,52,1311,53,1608,54,1608,55,1608,56,1608,57,1986,58,1986,59,1986,60,1986,61,2297,62,2297,63,2297,64,2297,65,2398,66,2398,67,2398,68,2398,69,2430,70,2430,71,2430,72,2430,73,2441,74,2441,75,2441,76,2441,77,2445,78,2445,79,2445,80,2445,81,2445,82,2445,83,2445]],
2446]
threadcount = len(data)-1
totaltime = data[-1]

fig, ax = plt.subplots()

threadcount = len(data)-1;
for i in range(0,threadcount):
    row = data[i][1]
    rowlen = len(row)
    lm = list()
    lq = list()
    for j in range(0,rowlen,4):
        piecenumber = row[j]//4
        code = row[j]%4
        start = row[j+1]
        stop = row[j+3]
        if code == 2:
            lq.append((start, stop - start))
            ax.text((start+stop)/2.0, i - 0.1, str(piecenumber),ha='center', va='center')
        else:
            lm.append((start, stop - start))
            ax.text((start+stop)/2.0, i + 0.3, str(piecenumber),ha='center', va='center')
    ax.broken_barh(lq, (i+0.0, 0.2), facecolors='grey')
    ax.broken_barh(lm, (i+0.4, 0.2), facecolors='lightgrey')

ax.set_ylim(-1, threadcount+1)
ax.set_xlim(0, totaltime)
ax.set_xlabel('ms since start (total time ' + str(totaltime) + 'ms)')
ax.set_yticks(list(range(0, threadcount)))
ax.set_yticklabels(list(range(0, threadcount)))
ax.grid(True)
ax.set_aspect((0.2*totaltime)/(threadcount+2));
plt.show()
