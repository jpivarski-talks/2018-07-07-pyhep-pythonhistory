Directory[]

Out[1]= /A/lns101/nfs/homes/cleo/mccann/public_html/private/my_gammaee/tools

<< "../fits/y1s_kf_dumpsave.mx"

<< "kf_fits.m"

MultipleListPlot[
    Table[ { { i, GetParams[y1sfits[[i,3]]][[1]] },
	     ErrorBar[Sqrt[GetErrorMatrix[y1sfits[[i,3]]][[1,1]]]] },
	   { i, 1, Length[y1sfits] } ]
    , Axes -> False
    , Frame -> True
    , FrameLabel -> { "Scan Number", "Area in MeV nb",
		      "Area Measurement for Each Scan", None }
    , ImageSize -> 72 * 5
		]
MultipleListPlot[
    Table[ { { i, GetParams[y1sfits[[i,3]]][[2]] },
	     ErrorBar[Sqrt[GetErrorMatrix[y1sfits[[i,3]]][[2,2]]]] },
	   { i, 1, Length[y1sfits] } ]
    , Axes -> False
    , Frame -> True
    , FrameLabel -> { "Scan Number", "Mass in GeV",
		      "Mass Measurement for Each Scan", None }
    , ImageSize -> 72 * 5
		]
MultipleListPlot[
    Table[ { { i, GetParams[y1sfits[[i,3]]][[3]] },
	     ErrorBar[Sqrt[GetErrorMatrix[y1sfits[[i,3]]][[3,3]]]] },
	   { i, 1, Length[y1sfits] } ]
    , Axes -> False
    , Frame -> True
    , FrameLabel -> { "Scan Number", "Width in MeV",
		      "Width Measurement for Each Scan", None }
    , ImageSize -> 72 * 5
		]
MultipleListPlot[
    Table[ { { i, GetParams[y1sfits[[i,3]]][[4]] },
	     ErrorBar[Sqrt[GetErrorMatrix[y1sfits[[i,3]]][[4,4]]]] },
	   { i, 1, Length[y1sfits] } ]
    , Axes -> False
    , Frame -> True
    , FrameLabel -> { "Scan Number", "Background in nb",
		      "Background Measurement for Each Scan", None }
    , ImageSize -> 72 * 5
		]

meanmean = GetParams[y1sidealfit][[2]]

Out[5]= 9.46007

shifts = Table[ GetParams[y1sfits[[i,3]]][[2]],
		{ i, 1, Length[y1sfits] } ] - meanmean

Out[6]= {-0.0000162542, -0.000291023, 0.000215811, 0.000168069, 0.0000526342, 
 
>    -0.000172138, -0.0000654832}

y1sDataShifted = Table[ { #[[1]] - shifts[[i]], #[[2]], #[[3]] }&
			/@ y1sData[[i]],
			{ i, 1, Length[y1sfits] } ];

PlotFuncData[ y1sidealfit, Flatten[ y1sDataShifted, 1 ],
	      "All Upsilon-1s Data to Kuraev-Fadin Fit",
	      { ImageSize -> 72 * 7 } ];

messofdata = Join[ Flatten[ y1sDataShifted, 1 ], Flatten[ y1sExtra, 1 ] ];

{ messchi2, messdof, messfunc } = FitToFunc[ y1sidealfit, messofdata ]
[Mathematica exited abnormally with code 1.]

[Mathematica segmentation fault.]

Exit
[Mathematica finished.]
