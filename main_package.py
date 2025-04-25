import numpy as np
import matplotlib.pyplot as plt

def main():
	return

class position_distribution:
    def __init__(self, file_pwd):
        self.data = np.loadtxt(file_pwd)
        self.dis_r = None
        self.number_molecules = len(self.data[0])
        self.BigO = [self.data[:, i:i+3] for i in range(0,self.number_molecules,3)]
    
    def sample_dis(self):
        distance = 25.2
        self.dis_r = []

        for j in range(int(self.number_molecules/3)):
            for k in range(j+1, int(self.number_molecules/3)):
                for i in range(np.shape(self.BigO)[1]):
                    x = np.abs(self.BigO[j][i][0] - self.BigO[k][i][0])
                    x = np.mod(x, distance)

                    if abs(x-distance) < x:
                        x = abs(x-distance)

                    y = np.abs(self.BigO[j][i][1] - self.BigO[k][i][1])
                    y = np.mod(y, distance)

                    if abs(y-distance) < y:
                        y = abs(y-distance)

                    z = np.abs(self.BigO[j][i][2] - self.BigO[k][i][2])
                    z = np.mod(z, distance)

                    if abs(z-distance) < z:
                        z = abs(z-distance)
                
                    r = np.sqrt(x**2+y**2+z**2)
                    
                    # if r <= rcut:
                    self.dis_r.append(r)

    def plot(self, bins): 
        plt.hist(np.array(self.dis_r), bins=bins, color = 'gray', edgecolor = 'black', density=True, alpha=0.5, label=r"O-O")
        plt.legend()
            

if __name__ == '__main__':
	main()
