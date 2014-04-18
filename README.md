Model to assess the operation possibilities of CHP's in cooperation with wind turbines
======================================================================================

Motivation
----------
I wrote this model for my master thesis in 2013-2014. 

Goal
----
The goal of the model is to reduce the profit loss in providing wind energy due to the imbalance prices. Therefore, residential CHP's are used. The model assumes N different houses that all have a certain heat demand and a CHP. The heat demand must be met at all times.

Using forecasts of wind production, imbalance prices, heat demand and the spotprice, an optimal bidding is done on the day-ahead market. Next, these results are compared with the most optimal results (the actuals - wind production, imbalance prices, heat demand and spotprice - are known) to assess the quality of the bidding.

In the model, parameters such as the wind turbine installed power, the number of units, the amount of generated scenarios and forecasts can be adapted.

Assumptions
-----------
* The model uses data from the Belgian TSO, Elia.

Howto
-----
1. The model is a combination of Matlab code (.m) and GAMS using the Cplex solver (.gms). These programs and packages are required to run the model.
2. Open the Matlab code: CHP.m
3. Play around with the variables under DEFINING VARIABLES
4. Run the script with Matlab
5. The results, together with the input parameters are stored in results.mat (different runs with different parameters are added to results.mat so previous results are not thrown away)
6. Stability.m shows a summary of all the results as dots on a graph

If you want to use any of the code, please send an email to jef(dot)daniels1(at)gmail(dot)com

Enjoy!
