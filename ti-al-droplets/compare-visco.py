import matplotlib.pyplot as plt

alpha = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

NVT_avg = [0.31950506989270705, 0.7160339798312625, 0.623996379994416, 0.5907776243106834, 1.3858058006604548, 1.2034065142697432]
NVT_err = [0.13897683837650715, 0.30819235284775837, 0.3257442951398426, 0.2634096889709532, 0.7308669264976674, 0.4894864243865616]

NPT_avg = [0.41830572279582734, 0.5518723973989429, 0.4617192498081739, 1.6056277782591855, 1.4317475503399009, 1.491262383295906]
NPT_err = [0.17956090534915953, 0.2459351737264289, 0.16956693620742622, 0.8660591734379591, 0.4330199992250701, 0.1938842554243243]

NVE_avg = [0.378634671712897, 0.333588225124117, 0.5589408538837858, 0.9035879477349777, 1.8292988391334315, 1.8188915080393622]
NVE_err = [0.24221355042670803, 0.26427286955633134, 0.2242471297337234, 0.1827315645014213, 0.8162510952031284, 0.8116072419288871]

eb1 = plt.errorbar(alpha, NVT_avg, yerr=NVT_err, fmt='ro', label="NVT", 
    ms=12, markerfacecolor='w', markeredgewidth=2.5, capsize=5)
eb1[-1][0].set_linestyle(':')
eb2 = plt.errorbar(alpha, NPT_avg, yerr=NPT_err, fmt='bs', label="NPT",
    ms=10, markerfacecolor='w', markeredgewidth=2.5, capsize=5)
eb2[-1][0].set_linestyle(':')
eb3 = plt.errorbar(alpha, NVE_avg, yerr=NVE_err, fmt='kD', label="NVE",
    ms=9, markerfacecolor='w', markeredgewidth=2.5, capsize=5)
eb3[-1][0].set_linestyle(':')

plt.legend(fontsize=20)
plt.ylabel(r"$\eta$ [cP]", fontsize=25)
plt.xlabel("Ti at%", fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.show()