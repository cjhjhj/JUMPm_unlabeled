#!/usr/bin/python

class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self):
        self.count += 1
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\rProgress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                   int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()

# class progressBar():
#     def __init__(self, nTot):
#         # Initialization
#         self.nTot = nTot
#         self.minorTicks = 5
#         self.majorTicks = 25
#         self.cycleInterval = 50
#         self.cycleCharArray = ["-", "/", '\\', '|']
#         self.cycleI = 0
#         self.curPercent = 0
#         self.counter = 0
#         self.backupCursor = 0
#
#     def increment(self):
#         self.counter += 1
#         if self.counter % self.cycleInterval == 0:
#             if self.backupCursor == 1:
#                 print("\b", end='')
#             print(self.cycleCharArray[self.cycleI], end='')
#             self.cycleI = (self.cycleI + 1) % len(self.cycleCharArray)
#             self.backupCursor = 1
#         if (self.counter / self.nTot * 100) >= self.curPercent:
#             if self.curPercent % self.majorTicks == 0:
#                 if self.backupCursor == 1:
#                     print("\b", end='')
#                 print("%d%%" % self.curPercent, end='')
#                 self.backupCursor = 0
#                 if self.curPercent == 100:
#                     print()
#             else:
#                 if self.backupCursor == 1:
#                     print("\b", end='')
#                 print(". ", end='')
#             self.curPercent += self.minorTicks
#
#
# nTot = 400;
# progress = progressBar(nTot)
# for i in range(0, nTot):
#     progress.increment()