### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ ee98f90e-7c5a-11eb-31bf-c396e378e6d9
using DataFrames, Dates, GeoJSON, GeoInterface, Glob, HTTP, JSON, NCDatasets, PlutoUI, PolygonOps, StatsPlots, Statistics

# ╔═╡ 04e64e8c-770d-11eb-200c-6b46ee645dc1
html"<button onclick=present()>Present</button>"

# ╔═╡ 6780f078-76dc-11eb-2e32-ed327bd5859d
md"# Case Study: Impact of the 2019 Midwest Flood on SIF over Illinois"

# ╔═╡ 795da832-76de-11eb-1731-078f75664e89
md"The US Midwest experienced a major flooding event in Spring 2019. Here, we will use TROPOMI and OCO-2 SIF to investigate how the seasonal cycle of photosynthesis has changed in response to the delayed planting. 

Further reading:

[NY Times article](https://www.nytimes.com/interactive/2019/09/11/us/midwest-flooding.html)

[Caltech article](https://www.caltech.edu/about/news/flooding-stunted-2019-cropland-growing-season-resulting-more-atmospheric-cosub2sub)

[Yin et al. 2020](https://doi.org/10.1029/2019AV000140)

[Source](https://eric.clst.org/tech/usgeojson/) of the GeoJSON file of US states
"

# ╔═╡ 3bb15b3e-7703-11eb-2395-3760c8750a3a
md"
This notebook uses the functions from the first one with minor modifications.

> We follow the same procedure to load the data from files as in the first notebook, so it will take a moment to loop through all the files."

# ╔═╡ d716d466-770c-11eb-01c4-dda250d53312
md"# Load required packages"

# ╔═╡ beba0132-76f9-11eb-09c4-0f2c50b96500


# ╔═╡ cfe90be2-76f9-11eb-0065-c94a4adb1303
plotly()

# ╔═╡ 08830152-7705-11eb-3032-51b754b5f526
md"# Reading the boundary coordinates"

# ╔═╡ 82b7c122-7705-11eb-039b-d1d99194940c
md"Reading the shapefile containing all US states:"

# ╔═╡ 9f2e54a4-7711-11eb-2016-fbf66af577f0
US_states=GeoJSON.read(String(HTTP.get("https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json").body))

# ╔═╡ 38a3d29a-76ff-11eb-19d3-33482fe4000d
#length(US_states.features)

# ╔═╡ dcc9d15c-76fe-11eb-398a-d5d431aee40e
md"Read the shape of Illinois:"

# ╔═╡ 2bbb67a4-76fe-11eb-2b30-4f7e4dbaab97
begin
	stateNames = []
	for i in 1:length(US_states.features)
		push!(stateNames,US_states.features[i].properties["NAME"])
	end
	#US_states.features[2].properties["NAME"]
	#idxIllinois = 
	Illinois = US_states.features[findfirst(stateNames .== "Illinois")[1]].geometry
end

# ╔═╡ 6c7f69ec-76ff-11eb-0aa9-2767374c75ba
plot(Illinois)

# ╔═╡ b1175210-7c89-11eb-15dd-afdb7d469af7
IllinoisLonLat = Tuple.(Illinois.coordinates[1])

# ╔═╡ fe40fcf4-7706-11eb-02c9-e9f8e2e1d59e
md"# Reading Data for Illinois"

# ╔═╡ 6f558260-7c5d-11eb-2f33-c7b983f49495
md" **Define sensor here:**"

# ╔═╡ 7c78f170-7c5d-11eb-3547-256ee2e742ea
sensor = "TROPOMI"
#sensor = "OCO-2"

# ╔═╡ 3f8b6092-7c5b-11eb-39bc-f7e3291e6aac
begin
	getFiles = function(sensor="TROPOMI") 
		## TROPOMI:
		if sensor == "TROPOMI"
			sifDataDir = "/net/fluo/data2/data/TROPOMI_SIF740nm/original"
			allFiles = glob("*.nc",sifDataDir)

			FilePattern = "TROPO_SIF_YYYY-MM-DD_*.nc"	
			posY = findfirst("YYYY",FilePattern)
			posM = findfirst("MM",FilePattern)
			posD = findfirst("DD",FilePattern)
			allDates=Vector{Date}()
			for file in basename.(allFiles)
				push!(allDates,Date(file[posY]*file[posM]*file[posD], "yyyymmdd"))
			end
			
			
		## OCO-2:
		elseif sensor == "OCO-2"
			sifDataDir = "/net/fluo/data1/ftp/data/OCO2/sif_lite_B8100"
			allFiles = glob("*/*/*.nc4",sifDataDir)
		
			FilePattern = "oco2_LtSIF_YYMMDD*.nc4"
			posY = findfirst("YY",FilePattern)
			posM = findfirst("MM",FilePattern)
			posD = findfirst("DD",FilePattern)
	
			allDates=Vector{Date}()
			for file in basename.(allFiles)
				push!(allDates,Date("20"*file[posY]*file[posM]*file[posD], "yyyymmdd"))
			end
		end
		
		return(allFiles, allDates)

	end
	
	
	getFileIdx4Year = function(yyyy)	
		startDate = Date(yyyy,1,1)
		stopDate  = Date(yyyy,12,31)
		FileInd = findall(x -> x >= startDate && x <= stopDate, allDates)
		return(FileInd)
	end
	
	
end

# ╔═╡ 251b84d0-7c5b-11eb-01eb-93007209c2e8
begin
	## Edited dictionary to load only necessary fields
	
	getBaseDict = function(sensor="TROPOMI")
		
		if sensor == "TROPOMI"
			
		baseDict = JSON.Parser.parse("""
    	    [{
				"basics": {
        		"lat": "lat",
        		"lon": "lon",
				"sif": "sif",
        		"time": "TIME"
				},
				"time_struc": "UTC, seconds since 1970-01-01 00:00:00"
			}]
		""")
			
		elseif sensor == "OCO-2"

		baseDict = JSON.Parser.parse("""
  	      [{
				"basics": {
  	 	    	"lat": "latitude",
      		  	"lon": "longitude",
				"sif": "SIF_757nm",
				"measurement_mode": "measurement_mode",
       		 	"time": "time",
				},
				"time_struc": "seconds since 1993-1-1 0:0:0"
			}]
		""")
		end
		
		return(baseDict[1]["basics"])
	end
	
	
	
	
end

# ╔═╡ d651f31a-7c5b-11eb-3432-21404096f38f
begin
	
numberNames = function(what::String,nWhat::Int)
	numNames = Vector()
	for x in string.(Vector(1:nWhat))
		push!(numNames,what*x)
	end
	return(numNames)
end

# now return DataFrame of filtered data:
readNCdata = function(NCFiles,baseDict,polyLonLat)
	
	## strategy: 
	##			- collect data for a minimum boundary box 
	##			- filter for soundings within polyLonLat and return a data frame	
	
	lons = [x[1] for x in polyLonLat]
	lats = [x[2] for x in polyLonLat]
	
	##mbr format: lonmin,lonmax,latmin,latmax
	mbr = [minimum(lons),maximum(lons),minimum(lats),maximum(lats)]
	
	#idxKeys = findall(in(["lon","lat"]).(string.(keys(baseDict))) .== false)	
	
	dfOut = DataFrame()
	
	for aFile in NCFiles
	
		dfTmp = DataFrame()
		## check first if the file has any soundings inside the mbr before reading the other variables (? -- TBD)
		#dfTmp.lon = Dataset(aFile)[baseDict["lon"]].var[:]
		#dfTmp.lat = Dataset(aFile)[baseDict["lat"]].var[:]	
		#for aKey in string.(keys(baseDict))[idxKeys]
		
		for aKey in keys(baseDict)
			tmp = Dataset(aFile)[baseDict[aKey]].var[:]
			if length(size(tmp))==1
				dfTmp[!,aKey] = tmp
			else
				## check which dimension is sounding and lat/lon (different for TROPOMI and OCO-2!)	
				if size(tmp)[2] > size(tmp)[1] 
					extendedKeys=numberNames(aKey,size(tmp)[1])
					for iKey in 1:size(tmp)[1]
						dfTmp[!,extendedKeys[iKey]] = tmp[iKey,:]
					end				
				else
					extendedKeys=numberNames(aKey,size(tmp)[2])
					for iKey in 1:size(tmp)[2]
						dfTmp[!,extendedKeys[iKey]] = tmp[:,iKey]
					end
				end
			end
		end
	
		# filtering:
		filter!("lon" => >=(mbr[1]), dfTmp)
		filter!("lon" => <=(mbr[2]), dfTmp)
		filter!("lat" => >=(mbr[3]), dfTmp)
		filter!("lat" => <=(mbr[4]), dfTmp)		

		if size(dfTmp)!=(0,0)
			points = []
			for i in 1:length(df.lon)
				push!(points, (df.lons[i],df.lats[i]) )
			end

			dfTmp.inside = [inpolygon(p, polyLonLat; in=true, on=true, out=false) for p in points]
			
			filter!("inside" => .==(true), dfTmp)
			## drop	"inside" column:
			select!(dfTmp, Not("inside"))
		end
			
		# attach to output data frame:
		if size(dfTmp)!=(0,0)
			if size(dfOut)==(0,0) 
				dfOut = dfTmp
			else
				append!(dfOut, dfTmp)
			end
		end
	
	end

	return(dfOut)

end
	
end

# ╔═╡ 2827c5b4-7c5c-11eb-1886-e323231c7c07
begin
	getOneYear = function(yyyy, polyLonLat, sensor="TROPOMI")
		
		allFiles, allDates = getFiles(sensor)
		baseDict = getBaseDict(sensor)
		FileInd  = getFileIdx4Year(yyyy)
		
		df = readNCdata(allFiles[FileInd],baseDict,polyLonLat)
		
		############# Convert time strings to julia's date type
		
		# TROPOMI:
		if sensor == "TROPOMI"
			df.utc = unix2datetime.(df.time)
		# OCO-2:
		elseif sensor == "OCO-2"
			filter!("measurement_mode" => ==(0), df)
			timeOffset = datetime2unix(DateTime(1993,1,1, 0,0,0))
			df.utc = unix2datetime.(df.time .+ timeOffset)
		end
		
		# Add also a column with plain days:
		df.Date = Date.(df.utc)
				
		return(df)
	end
	
end

# ╔═╡ 72bc1716-770d-11eb-3f50-e91939e494b9
begin
		#Illinois2018 = getOneYear(2018, polyLonLat=IllinoisLonLat, sensor=sensor)
		#Illinois2019 = getOneYear(2019, polyLonLat=IllinoisLonLat, sensor=sensor)
end

# ╔═╡ 9ea08f84-770c-11eb-31c7-850f3244c4ff
md"# Temporal Averaging"

# ╔═╡ bed15e68-770d-11eb-16ad-19c4dc89e4d8
md"> Now we see that the 2019 flood had a significant impact on the seasonality of photosynthesis over Illinois!"

# ╔═╡ 72dd46f6-7709-11eb-11e3-b7cccdda854b
begin
	time_slider= @bind n_weeks Slider(1:1:8; default=2)
	md"""
	**Averaging period (weeks):** $(time_slider)
	"""
end

# ╔═╡ 8c89d290-7709-11eb-2537-2db49c3d487b
md"**$(n_weeks) Week(s)** "

# ╔═╡ 46e18090-7708-11eb-0d74-f77ec2e38bc2
begin
#df.Date = Date.(df.utc)
	
## function returns temporal average (mean, std, n)	
tempAv = function(df::DataFrame, nWeeks, what::String)
	
	# rename the column of interest
	subDF = select(df, what => "what", "Date" => "Date")	
	
	# generate breaks:	
	tBreaks = collect(minimum(subDF.Date):Week(nWeeks):maximum(subDF.Date)+Day(1))
	
	## initialize avTime column:
	subDF.avTime = Array{DateTime, 1}(undef,  size(subDF)[1])
	
	## populate avTime:
	for i in 1:(length(tBreaks)-1)
			subDF.avTime[findall((subDF.Date .>= tBreaks[i]) .& (subDF.Date .< tBreaks[i+1]))] .= DateTime(tBreaks[i]) .+ (DateTime(tBreaks[i+1]) - DateTime(tBreaks[i])) ./ 2
	end
	## clean-up incomplete time steps:
	filter!("avTime" => >(minimum(tBreaks)), subDF)
	filter!("avTime" => <(maximum(tBreaks)), subDF)

	## generate output	
	outDF = combine(groupby(subDF, :avTime), :what => mean)
	rename!(outDF, :what_mean => :mean)
	outDF.sd = combine(groupby(subDF, :avTime), :what => std)[:,2]
	outDF.n   = combine(groupby(subDF, :avTime), :what => length)[:,2]
	#sort!(outDF, :avTime)
	return(outDF)	
end
	
end

# ╔═╡ 0aca52ac-7709-11eb-0414-311256298017
begin
	av2018 = tempAv(Illinois2018, n_weeks, "sif")	
	av2019 = tempAv(Illinois2019, n_weeks, "sif")	

	plot(Dates.dayofyear.(av2018.avTime),av2018.mean; ribbon=av2018.sd, linewidth=4, label= "2018 mean and sd ("*string(n_weeks)*" weeks)")
	plot!(Dates.dayofyear.(av2019.avTime),av2019.mean; ribbon=av2019.sd, linewidth=4, label= "2019 mean and sd ("*string(n_weeks)*" weeks)")
end

# ╔═╡ df940696-7c7b-11eb-0668-6df03990aab1
http://localhost:2345/?secret=qespzx0O

# ╔═╡ d55bdff0-7c7b-11eb-2017-6f417d68079e


# ╔═╡ Cell order:
# ╟─04e64e8c-770d-11eb-200c-6b46ee645dc1
# ╟─6780f078-76dc-11eb-2e32-ed327bd5859d
# ╟─795da832-76de-11eb-1731-078f75664e89
# ╟─3bb15b3e-7703-11eb-2395-3760c8750a3a
# ╟─d716d466-770c-11eb-01c4-dda250d53312
# ╟─beba0132-76f9-11eb-09c4-0f2c50b96500
# ╠═ee98f90e-7c5a-11eb-31bf-c396e378e6d9
# ╠═cfe90be2-76f9-11eb-0065-c94a4adb1303
# ╟─08830152-7705-11eb-3032-51b754b5f526
# ╟─82b7c122-7705-11eb-039b-d1d99194940c
# ╠═9f2e54a4-7711-11eb-2016-fbf66af577f0
# ╠═38a3d29a-76ff-11eb-19d3-33482fe4000d
# ╟─dcc9d15c-76fe-11eb-398a-d5d431aee40e
# ╠═2bbb67a4-76fe-11eb-2b30-4f7e4dbaab97
# ╠═6c7f69ec-76ff-11eb-0aa9-2767374c75ba
# ╠═b1175210-7c89-11eb-15dd-afdb7d469af7
# ╟─fe40fcf4-7706-11eb-02c9-e9f8e2e1d59e
# ╟─6f558260-7c5d-11eb-2f33-c7b983f49495
# ╠═7c78f170-7c5d-11eb-3547-256ee2e742ea
# ╠═3f8b6092-7c5b-11eb-39bc-f7e3291e6aac
# ╠═251b84d0-7c5b-11eb-01eb-93007209c2e8
# ╠═d651f31a-7c5b-11eb-3432-21404096f38f
# ╠═2827c5b4-7c5c-11eb-1886-e323231c7c07
# ╠═72bc1716-770d-11eb-3f50-e91939e494b9
# ╟─9ea08f84-770c-11eb-31c7-850f3244c4ff
# ╟─bed15e68-770d-11eb-16ad-19c4dc89e4d8
# ╟─8c89d290-7709-11eb-2537-2db49c3d487b
# ╟─72dd46f6-7709-11eb-11e3-b7cccdda854b
# ╠═0aca52ac-7709-11eb-0414-311256298017
# ╠═46e18090-7708-11eb-0d74-f77ec2e38bc2
# ╠═df940696-7c7b-11eb-0668-6df03990aab1
# ╟─d55bdff0-7c7b-11eb-2017-6f417d68079e
