----Tick-Host Model----

mainTick.m - run this file to find fixed points and 2-cycles of the tick-host system and determine local stability. This file generates Figure 1 our paper, which depicts the larval tick population with respect to the demographic reproductive number.

R0.m - run this file to calculate the basic reproductive number R_0 on the fixed points and 2-cycles. This file generates Figure 2 in our paper, which depicts R_0 with respect to the demographic reproductive number.

makeJacobianTick.m - this is called by the "mainTick.m" and "R0.m" files to assess local stability; it numerically finds the Jacobian of the tick-host system without infection.

makeJacobianInfection.m - this is called by the "R0.m" file to assess local stability; it numerically finds the Jacobian of the tick-host system with infection.

mySystemTick.m - this is called by the "mainTick.m" and "R0.m" files during the system solving process; it iterates the tickMap function once to find fixed points.

mySystemTick2.m - this is called by the "mainTick.m" and "R0.m" files during the system solving process; it iterates the tickMap function twice to find 2-cycle points.

tickMap.m - this is called by the "mySystemTick.m" and "mySystemTick2.m" files during the system solving process; describes the tick-host model without infection.

----Host-Parasitoid Model----

mainParasitoid.m - run this file to find fixed points, 2-cycles, 4-cycles, and 8-cycles of the host-parasitoid system and determine local stability. This file generates Figure 3 our paper, which depicts the host and parasitoid populations with respect to the demographic reproductive number.

makeJacobianParasitoid.m - this is called by the "mainParasitoid.m" file to assess local stability; it numerically finds the Jacobian of the host-parasitoid system without infection.

mySystemParasitoid.m - this is called by the "mainParasitoid.m" file during the system solving process; it iterates the parasitoidMap function once to find fixed points.

mySystemParasitoid2.m - this is called by the "mainParasitoid.m" file during the system solving process; it iterates the parasitoidMap function twice to find 2-cycle points.

mySystemParasitoid4.m - this is called by the "mainParasitoid.m" file during the system solving process; it iterates the parasitoidMap function four times to find 4-cycle points.

mySystemParasitoid8.m - this is called by the "mainParasitoid.m" file during the system solving process; it iterates the parasitoidMap function eight times to find 8-cycle points.

parasitoidMap.m - this is called by the "mySystemParasitoid.m", "mySystemParasitoid2.m", "mySystemParasitoid4.m" and "mySystemParasitoid8.m" files during the system solving process; describes the host-parasitoid model without infection.
