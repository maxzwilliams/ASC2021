In this folder all the codes are stored.

WARNINGS (DO NOT IGNORE):
LB.go generates many large .csv files that can easily fill up your computers storage. If left unattended, LB.go can fill up your computers storage. THIS IS VERY BAD, YOU MAY NOT BE ABLE TO LOG BACK INTO YOUR COMPUTER (this happened to me)
Generating plots from c.py, p.py and plots.py can and will fill up your computers memory after a short time (a few minutes for 8GiB of memory). This will eventually crash your computer, but wont screw it up after you reboot or it crashes.


There are four codes, p.py, c.py,  go/LB.go, go/plot.py

p.py is the code that simulates the polar coordinate streamfunction code
c.py is the code that simulates the cartesian coordinate streamfunction code

go/LB.go is the code that works by Lattice Boltzman
go/plot.py plots data produced by LB.go


How to run p.py and c.py:
To run p.py or c.py navigate to the codes folder in your terminal. Then run "python3 p.py" or "python3 c.py". The program should run if you have python3 installed. Note that you may need to use py or python or py3 inplace of python3

While running c.py will write plots of temperature, vorticity and streamfunction in the folders t, vts, sfs. Similary, p.py will write temperature, vorticity and streamfunction to tpolar, vtspolar, sfspolar.

How to run LB.go and render its outputs:
First navigate to the go folder in the terminal. Run "go run LB.go" and the program LB.go should run. While running, it will write csv files that detail the internal energy, density and velocity fields into the folders ep, rho and u. 

While LB.go is running, or once its done (note this will take a very long time at default settings > 1 week) you can render the data into plots using plot.py.
Open a terminal and navigate to go. Now run "python3 plot.py <startIndex> <endIndex> <stepSize>". This will render plots of timestep <startIndex> to <endIndex> with a stepsize of <stepSize>. Note that LB.go doesnt output data at every timestep and so a stepSize that is too small will either fail to find the required file or may read some other file from an older running of the program. LB.go records every 100th iteration. So I would recomend running "python3 plot.py 0 <endIndex> 100" where <endIndex> is some multiple of 100. Once plot.py is complete it will write plots of the internal energy, density and velocity fields into pep, prho, pu.





