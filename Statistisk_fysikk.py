import numpy as np
from numpy import linalg as LA

import matplotlib.pyplot as plt


#Generate alternating list of ones and zeros

nMax = 4

n = 5

B_ex = 1 #vet ikke hvor høy den skal være.

J_per, J_par , beta_const = 1, 1, 0.5

Beta_range = np.linspace(0.0005, 1.0, 100)


sigma_matrix = np.ones([2**n, n])

for i in range(n):
    sigma_matrix.T[i] = np.array([[1]*2**(i) + [-1]*2**(i)] * ((2**n)//(2**(i+1)))).flatten()

#print(sigma_matrix)



def P_matrix(beta_const):

    P = np.ones([2**n, 2**n])

    for l in range(2**n):
        for m in range(2**n):
            condition = l
            if 2**n == l + 1:
                condition = 0
                P[l][m] = np.exp(beta_const * B_ex * J_par * sum(sigma_matrix[l] * sigma_matrix[m] + J_per*beta_const*(sigma_matrix[l] + sigma_matrix[condition]) + ((B_ex*beta_const)/2) * (sigma_matrix[l] + sigma_matrix[m])))
            else:
                P[l][m] = np.exp(beta_const * B_ex * J_par * sum(
                    sigma_matrix[l] * sigma_matrix[m] + J_per * beta_const * (sigma_matrix[l] + sigma_matrix[condition]) + ((B_ex * beta_const) / 2) * (sigma_matrix[l] + sigma_matrix[m])))


    return P





eigenvalues = np.zeros((len(Beta_range), 2**n))
for b in range(len(Beta_range)):
    eigenvalues[b] = np.linalg.eigvals((P_matrix(Beta_range[b])))

print(np.real(eigenvalues))


eigenvalues = np.real(eigenvalues.T)

for i in range(len(eigenvalues)):
    plt.semilogy(Beta_range, eigenvalues[i])
plt.show()

tr = eigenvalues.T
list = []
for i in range(len(eigenvalues)):
    list = np.amax(eigenvalues[i])

N_list = np.log(list)
print(N_list)
lambdas = []

for i in range(len(N_list)):
    lambdas = (N_list[i + 1] - N_list[i])/(Beta_range[i + 1] - Beta_range[i])

for i in range(len(lambdas)):
    plt.plot(Beta_range, lambdas)
plt.show()
print(lambdas.shape)

