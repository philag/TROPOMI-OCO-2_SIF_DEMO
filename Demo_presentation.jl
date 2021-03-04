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

# ╔═╡ c25b5f6c-61a4-11eb-2c34-cbf860b5d314
using DataFrames, Dates, Dierckx, Distributions, GeoJSON, GeoInterface, Glob, HTTP, JSON, NCDatasets, PlutoUI, PolygonOps, StatsPlots, Statistics

# ╔═╡ 6659c6d4-7182-11eb-3b65-d939a40b4325
html"<button onclick=present()>Present</button>"

# ╔═╡ 1c2705ee-61a4-11eb-33f9-835cc5caa8a8
md"# TROPOMI/OCO-2 SIF DEMO"

# ╔═╡ b7921386-7186-11eb-20d2-d973f7940a35
md" **In this notebook we will:**
> - **Read L2 satellite data of solar-induced chlorophyll fluorescence (SIF) inferred from TROPOMI and OCO-2**

> - **Create time series over spatial features**

> - **Generate spatial composites**

> - **Evaluate the measurement uncertainties**"

# ╔═╡ e7799730-7187-11eb-02a8-e7f76ed6081c
md"## Start-up Julia"

# ╔═╡ b434a196-61a4-11eb-19fe-976e23dec488
md"**Loading required packages** which provide necessary functionality (e.g., reading certain data types)"

# ╔═╡ a4854224-7184-11eb-291a-e3667b3c33e2
md"**Switch to Plotly backend** to make plots interactive:"

# ╔═╡ 5250d5d2-718f-11eb-03ab-5d29ef4ce6ad
plotly()

# ╔═╡ 34be05d8-7185-11eb-3bd3-d1e24885ed75
md"## Download data"

# ╔═╡ 8cc541f6-7185-11eb-3034-675a8bf6ff7c
md"
> **TROPOMI data can be downloaded from our [ftp server](ftp://fluo.gps.caltech.edu/data/tropomi)** "

# ╔═╡ 2be873b4-718d-11eb-1a33-a742dd4610ff
LocalResource("./src/images/data_calendar.png", :width => 400)

# ╔═╡ bb5cb8a2-718d-11eb-3694-cbe3adf78a8d
md"
> **OCO-2 data can be downloaded from our [ftp server](ftp://fluo.gps.caltech.edu/data/OCO2/sif_lite_B8100/) or [GES DISC](https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_SIF_10r/summary?keywords=oco2%20sif%20lite)** "

# ╔═╡ a33e9ac2-718c-11eb-18d4-b90c627f52ec
md"## After Downloading the data"

# ╔═╡ 8f48bb84-718e-11eb-333e-8f1f2df863d9
md" Satellite data is often shared via .NC (netCDF, network Common Data Form) files. OCO-2/3 and TROPOMI data is stored in daily files, while TROPOMI files on our ftp server are bundled to monthly *.tar.gz archives. Unpacking the data is in this case the first step: ```tar -xvf *.tar.gz```
"

# ╔═╡ 38f5b8de-7190-11eb-0437-497aa85610f0
md"> **Satellite data sets have typically a large volume. In case of TROPOMI, one day includes about 4M observations over land with a size of up to 150MB (one month ~3.5GB).**"

# ╔═╡ 2df3dee4-61a4-11eb-2003-5b56cfb49d01
md"**Define location of the downloaded and unpacked files:**"

# ╔═╡ 5b5b38b4-61a4-11eb-284b-e329267e0fb2
## TROPOMI:
#sifDataDir = "/net/fluo/data2/data/TROPOMI_SIF740nm/original"
## OCO-2:
sifDataDir = "/net/fluo/data1/ftp/data/OCO2/sif_lite_B8100"

# ╔═╡ 4f00b9e0-61a4-11eb-0d33-f3ec82cd5396
md"**Get a list of available files:**"

# ╔═╡ 2f0d9a3a-61a5-11eb-2fc1-650e465800c8
## TROPOMI:
#allFiles = glob("*.nc",sifDataDir)
## OCO-2:
allFiles = glob("*/*/*.nc4",sifDataDir)

# ╔═╡ 62fa453c-61a5-11eb-0b39-015a9ff4155c
md" **Convert filenames to dates:**"

# ╔═╡ a18e5c16-61c3-11eb-37cf-edd1266f9dde
begin
	getTROPOMIdates = function(allFiles)
		## TROPOMI:
		FilePattern = "TROPO_SIF_YYYY-MM-DD_*.nc"	
		posY = findfirst("YYYY",FilePattern)
		posM = findfirst("MM",FilePattern)
		posD = findfirst("DD",FilePattern)
		allDates=Vector{Date}()
		for file in basename.(allFiles)
			push!(allDates,Date(file[posY]*file[posM]*file[posD], "yyyymmdd"))
		end
		return(allDates)
	end

	getOCO2dates = function(allFiles)
		FilePattern = "oco2_LtSIF_YYMMDD*.nc4"
		posY = findfirst("YY",FilePattern)
		posM = findfirst("MM",FilePattern)
		posD = findfirst("DD",FilePattern)
	
		allDates=Vector{Date}()
		for file in basename.(allFiles)
			push!(allDates,Date("20"*file[posY]*file[posM]*file[posD], "yyyymmdd"))
		end
		return(allDates)
	end


end

# ╔═╡ 179972c6-7257-11eb-3a3d-45c672612fb3
allDates = getOCO2dates(allFiles)

# ╔═╡ bc20f27e-61d1-11eb-2497-179d53ab2360
md" **Select date range:**"

# ╔═╡ dfaee354-61d1-11eb-2a88-8bb2002c2136
begin
	#startDate = Date(2020,1,1)
	#stopDate  = Date(2020,1,15)
	startDate = Date(2019,1,1)
	stopDate  = Date(2019,12,31)
	FileInd = findall(x -> x >= startDate && x <= stopDate, allDates)
end

# ╔═╡ a8dc819e-7188-11eb-2c1b-8998b48d99f5
md"## Create and read spatial features"

# ╔═╡ cfd6f918-7189-11eb-33f4-c5133381580b
md" **Defining Regions of interest:**

> One simple way to search for interesting features is using [Google Earth](https://www.google.com/earth/).
"

# ╔═╡ c475fc80-718a-11eb-2c79-3f45c09aaca4
LocalResource("./src/images/google_screenshot.png")

# ╔═╡ 2680b0a4-718c-11eb-2bad-137068539830
md"## Converting Google's KML files"

# ╔═╡ 471fe14a-718c-11eb-23a2-a3c5ed97e288
md" The Keyhole Markup Language (**KML**) uses a tag-based structure with nested elements and attributes consistent with the Extensible Markup Language (**XML**)."

# ╔═╡ a748c364-71a2-11eb-21fb-2ddc00692957
LocalResource("./src/images/kml_screenshot.png")

# ╔═╡ a8564f6a-71a2-11eb-1797-67abf06ac8c1
md"
Markup languages generally define a set of rules for encoding documents in a format that is both human-readable and machine-readable. Another widely used format with this property is the JavaScript Object Notation (**JSON**) and **GeoJSON**, an open standard for encoding geographic data structures."


# ╔═╡ d753c838-71a2-11eb-151d-6394549ab869
LocalResource("./src/images/json_screenshot.png")

# ╔═╡ d7f995ce-71a2-11eb-3e35-b113a803a9fb
md"
There are plenty of tools like this [one](https://mygeodata.cloud/converter/kml-to-geojson) to convert **KML** to **GeoJSON** files, which we will use in the following."

# ╔═╡ 91b73620-7193-11eb-070f-573411f9e611
md"## Reading GeoJSON files"

# ╔═╡ df1138ee-61d5-11eb-34b6-67ce049efb74
testROI=GeoJSON.read(read("./src/TROPO-DEMO.json"))

# ╔═╡ 14e5f60c-6687-11eb-37f6-fbf0d9061530
testROI.features#[1].properties["name"]

# ╔═╡ ad814292-65d9-11eb-1390-35f78905835c
### Function to get coordinates w/o altitude for plotting and to make use of "inpolygon" from PolygonOps.jl later on:
getLonLat = function(x::Feature)
	LonLat = []
	for i in 1:length(x.geometry.coordinates[1])
		LonLat = push!(LonLat,x.geometry.coordinates[1][i][1:2])
	end
	return(Tuple.(LonLat))
end

# ╔═╡ e083f3fc-65da-11eb-290c-8df8464a4ff1
begin
	vegPolygon = getLonLat(testROI.features[1])
	Plots.plot(vegPolygon, legend=false)
	refPolygon = getLonLat(testROI.features[2])
	Plots.plot!(refPolygon, legend=false)
	#Plots.plot!(mbr)
end

# ╔═╡ 53d0a158-62a6-11eb-3f8b-af81223169ec
md"## Compute the minimum boundary box (mbr)"

# ╔═╡ 213784d6-7198-11eb-16a2-254812fbed1d
md" > This step is necessary to narrow down the huge data volume."

# ╔═╡ b88cdc64-6509-11eb-24c7-fd78be2de33c
# accessing the coordinates works like this:
size(testROI.features[1].geometry.coordinates[1])

# ╔═╡ 80cc3a5c-629f-11eb-3807-f51439287448
getMBReps = function(x::FeatureCollection, epsLatLon=0.2)
	# MBR: minimum boundary box
	box   = zeros(4)
 	for ifeat in 1:length(x.features)
	
		tmp = x.features[ifeat].geometry.coordinates[1] |> DataFrame
		
		if ifeat==1
			box[1] = minimum(tmp[1,:])
			box[2] = maximum(tmp[1,:])
			box[3] = minimum(tmp[2,:])
			box[4] = maximum(tmp[2,:])
		end
		
		if ifeat > 1
			if minimum(tmp[1,:]) < box[1]
				box[1] = minimum(tmp[1,:])
			end
			if minimum(tmp[2,:]) < box[3] 
				box[3] = minimum(tmp[2,:])
			end
			if maximum(tmp[1,:]) > box[2] 
				box[2] = maximum(tmp[1,:])
			end
			if maximum(tmp[2,:]) > box[4] 
				box[4] = maximum(tmp[2,:])
			end
		end
	
	end
	box[1] = box[1]-epsLatLon
	box[2] = box[2]+epsLatLon
	box[3] = box[3]-epsLatLon
	box[4] = box[4]+epsLatLon

	return(box)
end

# ╔═╡ b1b2bb48-64fd-11eb-1ce8-274175a2d8ff
mbr=getMBReps(testROI)

# ╔═╡ b3fa469c-61d4-11eb-17ca-31a9319d9560
md"## Defining which parts of the NC files are relevant"

# ╔═╡ 50938100-62af-11eb-0027-2bc902261db2
# look at what's available in the NC files:
let
	ds = Dataset(allFiles[1])#["sif"].var[:]
end

# ╔═╡ 6288c4c6-6585-11eb-24c3-4b39125563bb
#timeStruc = ds["TIME"].attrib["units"]

# ╔═╡ 35ce9bee-62a7-11eb-0118-9528ab9f33d8
begin
	
	tropoDict = JSON.Parser.parse("""
        [{
			"basics": {
        	"lat": "lat",
        	"lon": "lon",
        	"lat_bnd": "lat_bnds",
        	"lon_bnd": "lon_bnds",
			"sif": "sif",
			"sif_err": "sif_err",
			"nir": "NIR",
        	"time": "TIME"
			},
			"time_struc": "UTC, seconds since 1970-01-01 00:00:00"
		}]
	""")

	oco2Dict = JSON.Parser.parse("""
        [{
			"basics": {
  	     	"lat": "latitude",
        	"lon": "longitude",
        	"lat_bnd": "footprint_vertex_latitude",
        	"lon_bnd": "footprint_vertex_longitude",
			"sif": "SIF_757nm",
			"sif_err": "SIF_757nm_uncert",
			"measurement_mode": "measurement_mode",
        	"time": "time",
			"nir": "continuum_radiance_757nm"
			},
			"time_struc": "seconds since 1993-1-1 0:0:0"
		}]
	""")
	
end

# ╔═╡ 254018d0-62b7-11eb-2b1d-138ad8d088b1
## TROPOMI:
#baseDict = tropoDict[1]["basics"]
## OCO-2:
baseDict = oco2Dict[1]["basics"]

# ╔═╡ 82ade7fe-62b7-11eb-3aa4-e9b8ba1f5a58
#tmp = Dataset(allFiles[1])[baseDict["lat"]].var[:]

# ╔═╡ 50f46e82-64c3-11eb-29d9-933218ab2e66
md"## Create a DataFrame from NC files"

# ╔═╡ b30830b2-71a3-11eb-2699-53e1f7bcdb43
md"
> [DataFrames.jl](https://dataframes.juliadata.org/stable/) makes it easier to work with tabular data. Its design and functionality are similar to those of pandas in Python and data.frame/data.table in R.
"

# ╔═╡ a33e67bc-650a-11eb-2ad2-43fe662d8a78
md"**Collecting data for the specified time period:**"

# ╔═╡ 00f8867e-64f3-11eb-18f0-2f31200f7628
#df = readNCdata(allFiles[FileInd],baseDict,mbr)

# ╔═╡ 7ff0299c-7265-11eb-326b-e5659d2f2f58
#filter!("measurement_mode" => ==(0), df)

# ╔═╡ 7563ebee-64f5-11eb-30ca-f5e9d7856dd7
begin
	
numberNames = function(what::String,nWhat::Int)
	numNames = Vector()
	for x in string.(Vector(1:nWhat))
		push!(numNames,what*x)
	end
	return(numNames)
end

# now return DataFrame of filtered data:
readNCdata = function(NCFiles,baseDict,mbr)
	
	dfOut = DataFrame()
	
	for aFile in NCFiles
	
		dfTmp = DataFrame()			
		for aKey in keys(baseDict)
			tmp = Dataset(aFile)[baseDict[aKey]].var[:]
			if length(size(tmp))==1
				dfTmp[!,aKey] = tmp
			else
				## check which dimension is sounding and lat/lon	
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
	# attach to final output data:
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

# ╔═╡ eea1155e-64f5-11eb-3eed-f34a6000d5ea
#test numberNames function:
#numberNames("lat",4)

# ╔═╡ a3b5d552-650e-11eb-0db6-cdd12b877af8
md"**Convert time structure:**"

# ╔═╡ c21d0ec8-71a4-11eb-1aa7-b1f03c83a3b7
md" For **TROPOMI** it is unix time, so the conversion to a DateTime object is simple:"

# ╔═╡ a12da5e2-71a4-11eb-3f91-437744720d0d
tropoDict[1]["time_struc"]

# ╔═╡ eb55f0c2-6583-11eb-2572-1b96530ce77a
#df.utc = unix2datetime.(df.time)

# ╔═╡ f62cce7c-725a-11eb-3b3b-ad594257cda6
md" For **OCO-2** the time structure has to be converted in a slightly different way:"

# ╔═╡ f33daa04-725a-11eb-2831-6fc756438afe
oco2Dict[1]["time_struc"]

# ╔═╡ 724200f0-7262-11eb-1b6a-3723b249389d
begin
	#timeOffset = datetime2unix(DateTime(1993,1,1, 0,0,0))
	#df.utc = unix2datetime.(df.time .+ timeOffset)
end

# ╔═╡ f3550094-71a5-11eb-1bd6-8150f1ad5604
md"Add also a column with plain days for convenience:"

# ╔═╡ b99cff74-650e-11eb-3373-1d52b025e9a5
#df.Date = Date.(df.utc)

# ╔═╡ 4a8de51e-65db-11eb-2e99-c70cb7193464
md"## What percentage of each footprint covers the spatial polygon?"

# ╔═╡ 26556eb4-664f-11eb-16f6-2b8c0c1a5945
## basic test for append!:
#begin
#	a = [1.]
#	#b = collect(1:10)*0.1
#	b=append!(a,(1. .+ collect(0:10) * 0.1))
#	c=append!(a,(2. .+ collect(0:10) * 0.5))
#	
#end

# ╔═╡ ac7223a8-668d-11eb-1906-a51e5b2ce8c1
md"A few required functions:"

# ╔═╡ 9555b186-664d-11eb-088a-bd3beb009eab
begin

# function to return formatted vertices:
getSdngVerts = function(df::DataFrame, i::Int64)
	# may need to edit this function later to order the vertices automatically	
	# e.g. in case LatLons are mixed up in the original data set
	vertLat = [df.lat_bnd1[i],df.lat_bnd2[i],df.lat_bnd3[i],df.lat_bnd4[i]]
	vertLon = [df.lon_bnd1[i],df.lon_bnd2[i],df.lon_bnd3[i],df.lon_bnd4[i]]
	return(DataFrame(lon=vertLon,lat=vertLat))
end
	
# To convert sounding to equally spaced spatial points we need two functions:
divLine = function(lon1,lat1,lon2,lat2,n::UInt8)
		dLat = (lat2-lat1)/(n-1)
		dLon = (lon2-lon1)/(n-1)
		
		lats = append!([lat1], (lat1 .+ collect(1:(n-1)) * dLat))
		lons = append!([lon1], (lon1 .+ collect(1:(n-1)) * dLon))
		return(lons, lats)
end
	
# Now divide each rectangle into n points:
getPoints = function(sdngVerts::DataFrame, n::UInt8)
		vertLon = sdngVerts.lon
		vertLat = sdngVerts.lat
		# Get reference points for two lines at the extremes:
	  	line1lon, line1lat = divLine(vertLon[1],vertLat[1],vertLon[2],vertLat[2],n)
	  	line2lon, line2lat = divLine(vertLon[4],vertLat[4],vertLon[3],vertLat[3],n)
	    # Now generate the inner points
		LonLat = []
		for i in 1:n
			tmpLon, tmpLat = divLine(line1lon[i], line1lat[i],  
						   			 line2lon[i], line2lat[i], n)
			for i in 1:length(tmpLon)[1]
				push!(LonLat,(tmpLon[i], tmpLat[i]))
			end
		end
		
		return(Tuple.(LonLat))
end

	
# functions to loop through soundings in DataFrame (df) to return coverage of polygon (percentage if n=10)
compCoverage = function(polygon::Array{Tuple{Float64,Float64},1}, sdngVerts::DataFrame, n::UInt8)
	testPoints = getPoints(sdngVerts, n)
	inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in testPoints]
	return(length(inside[inside .==true]))
end

getCoverage = function(df::DataFrame, polygon::Array{Tuple{Float64,Float64},1})
	[compCoverage(polygon, getSdngVerts(df,i), UInt8(10)) for i in 1:size(df)[1]]
end
	
	
end

# ╔═╡ c88c215e-6650-11eb-35fe-4dbe87cfb463
## Just to test the new functions:
# 
#getSdngVerts(df,1)
#
#divLine(Float32(0.),Float32(0.),Float32(10.), Float32(10.), UInt8(4))
#
#getPoints(getSdngVerts(df,1), UInt8(10))
#
#getCoverage(vegPolygon, getSdngVerts(df,25), UInt8(10))

# ╔═╡ a5b67f32-668f-11eb-090f-07b24fe0660d
md"**Attach coverage of all polygons to DataFrame:**"

# ╔═╡ 964d3d38-665d-11eb-1eb4-adbb37c289b6
begin
	# 
	df.VegCov = getCoverage(df,vegPolygon)
	df.RefCov = getCoverage(df,refPolygon)
end

# ╔═╡ 98684c54-6660-11eb-175f-c96664cc2f04
begin
# how does the coverage per area looks like?
#	Plots.plot(df.VegCov)
#	Plots.plot!(df.RefCov)
end

# ╔═╡ 0ff6e15e-6661-11eb-3e48-3b73b977b0fb
md"## Let's take a first look at the data:"

# ╔═╡ 0cda4bd4-6691-11eb-20c7-9b5097e1e7e9
begin
	minCovVeg_slider= @bind minCovVeg Slider(0:10:100; default=50)
	minCovRef_slider= @bind minCovRef Slider(0:10:100; default=50)
	md"""
	**Vegetation:** $(minCovVeg_slider)
	
	**Reference:**  $(minCovRef_slider)
	"""
end

# ╔═╡ 388e8398-6692-11eb-16c8-2bcccd929b61
md"**Soundings in time series have a coverage of at least $(minCovVeg)% (Vegetation) and $(minCovRef)% (Reference)**"

# ╔═╡ 9ba7d79a-6692-11eb-30e8-8d548156f71c
begin
	indVeg  = findall(df.VegCov .>= minCovVeg)
	ntotVeg = length(findall(df.VegCov .> 0))#size(df)[1] 

	Plots.scatter(df.utc[indVeg],df.sif[indVeg], markeralpha=0.5, label= "Vegetation (n="*string(length(indVeg))*"/"*string(ntotVeg)*")")

	indRef = findall(df.RefCov .>= minCovRef)
	ntotRef = length(findall(df.RefCov .> 0))#size(df)[1] 
	Plots.scatter!(df.utc[indRef],df.sif[indRef], markeralpha=0.5, label= "Reference (n="*string(length(indRef))*"/"*string(ntotRef)*")")
	
end

# ╔═╡ 46c0b6a0-6824-11eb-2407-a37dca8ef4bf
md"## Temporal averaging"

# ╔═╡ 2140d9a2-682c-11eb-0be0-8d4cba6681bb
begin
	time_slider= @bind n_days Slider(1:1:31; default=7)
	md"""
	**Averaging period (days):** $(time_slider)
	"""
end

# ╔═╡ 6988f26c-682c-11eb-3261-ed8a7c839fb4
md"**$(n_days) Day(s)** "

# ╔═╡ e28b9fc0-6a57-11eb-183f-6f47f1c39d6a
begin
#df.Date = Date.(df.utc)
	
## function returns temporal average (mean, std, n)	
tempAv = function(df::DataFrame, nDays, what::String)
	
	# rename the column of interest
	subDF = select(df, what => "what", "Date" => "Date")	
	
	# generate breaks:	
	tBreaks = collect(minimum(subDF.Date):Day(nDays):maximum(subDF.Date)+Day(1))
	
	## initialize avTime column:
	subDF.avTime = Array{DateTime, 1}(undef,  size(subDF)[1])
	
	## populate avTime:
	for i in 1:(length(tBreaks)-1)
			subDF.avTime[findall((subDF.Date .>= tBreaks[i]) .& (subDF.Date .< tBreaks[i+1]))] .= DateTime(tBreaks[i]) .+ (DateTime(tBreaks[i+1]) - DateTime(tBreaks[i])) ./ 2
	end
	## clean-up incomplete time steps:
	#filter!(row -> row[:avTime]!= DateTime("0000-12-31T00:00:00"), subDF)
	#filter!(row -> row[:avTime] > minimum(tBreaks), subDF)
	#filter!(row -> row[:avTime] < maximum(tBreaks), subDF)
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

# ╔═╡ 5623b812-6825-11eb-091d-3b676411207c
begin
av = tempAv(df[indVeg,:], n_days, "sif")	
scatter(df.utc[indVeg],df.sif[indVeg], markeralpha=0.2, c=:gray, label="single soundings" )

plot!(av.avTime,av.mean; ribbon=av.sd, linewidth=4, label= "mean and sd ("*string(n_days)*" days)")
end

# ╔═╡ d51c868a-760d-11eb-3420-9b51e5d9299a
av

# ╔═╡ e53f342c-669a-11eb-3afb-fbeede548218
md"## Computing spatial composites"

# ╔═╡ 59c35358-682d-11eb-0fd2-77744637cdbe
md"""**Resolution:**
		$(@bind resolution Select(["0.2" => "0.2°","0.1" => "0.1°","0.05" => "0.05°", "0.01" => "0.01°"]))
		"""

# ╔═╡ 2e2828e8-67ee-11eb-1b37-3ffc4a162aeb
begin
	
initGrid = function(resolution::Float64, mbr::Array{Float64,1}) ##mbr format: lonmin,lonmax,latmin,latmax
		
	lons = collect(mbr[1]+resolution/2.:resolution:mbr[2]-resolution/2.)
	lats = collect(mbr[3]+resolution/2.:resolution:mbr[4]-resolution/2.)
	
    #lons = collect(mbr[1]:resolution:mbr[2])
	#lats = collect(mbr[3]:resolution:mbr[4])
	
	iLons = collect(1:length(lons))
	iLats = collect(1:length(lats))
	
	## to map lon/lat to iLon/iLat
	iLon_ = Spline1D(lons,iLons)
	iLat_ = Spline1D(lats,iLats)

	grid = Array{Union{Missing,Float32}}(undef, length(lats), length(lons))
	grid[:,:] .= 0.
	
	return(lons, lats, iLon_, iLat_, grid)	
end

## generate polygon from corner coordinates:	
getPoly = function(xs,ys)
##first and last Tuple neeed to be identical to close the polygon
		polygon = []
		for i in 1:length(xs)
			push!(polygon, (xs[i], ys[i]))
		end
		push!(polygon, (xs[1], ys[1]))
	
		return(Tuple.(polygon))
end

	
oversample = function(df::DataFrame, resolution::Float64, mbr::Array{Float64,1}, what::String)
	
	#initialize grid:	
	lons, lats, iLon_,iLat_,grid = initGrid(resolution, mbr)	
	nGrid = copy(grid)
	
	## loop through soundings:	
	for iSdng in 1:size(df)[1]
		
		sdngVert = getSdngVerts(df,iSdng)
		
		# convert from coordinates to XY
		xs = iLon_(sdngVert.lon) 
		ys = iLat_(sdngVert.lat)	
		
		polygon = getPoly(xs,ys)
		
		# loop through all gridboxes in that area: 
		xmin = Int64(floor(minimum(xs)))
		if xmin==0 
			xmin += 1
		end
		xmax = Int64(ceil(maximum(xs)))
		if xmax > length(lons)
			xmax -= 1
		end
		ymin = Int64(floor(minimum(ys)))
		if ymin==0 
			ymin += 1
		end
		ymax = Int64(ceil(maximum(ys)))
		if ymax > length(lats)
			ymax -= 1
		end
			
		# compute coverage only if grid box vertices are not covered by sounding
		for x in xmin:xmax
			for y in ymin:ymax
				
				testPoints = [(x,y),(x,y+1),(x+1,y+1),(x+1,y)] # clockwise
				
				inside = [inpolygon(p, polygon; in=true, on=true, out=false) for p in testPoints]	
				## easy case if all corners are within sounding:
				if length(inside[inside .==true])==4
					nGrid[y,x] += 1.
					grid[y,x]  += df[iSdng,what]
				else	
				## compute coverage if gridcell is only partially covered:
					inflatedPoints = getPoints(DataFrame(lon=Float32.([x, x, x+1, x+1]), lat=Float32.([y, y+1, y+1, y])), UInt8(10))
					cov = [inpolygon(p, polygon; in=true, on=true, out=false) for p in inflatedPoints]
					weight = Float32(length(cov[cov .==true])/100.)
					nGrid[y,x] += weight
					grid[y,x]  += df[iSdng,what] * weight
					
				end
			end
		end
	end
	nGrid[nGrid .==0] .= missing
	avGrid = grid ./ nGrid

	return(avGrid, lons, lats)
end
	
	
end

# ╔═╡ 4a72d39c-65da-11eb-07c2-2dd1c2bcac1e
begin
	### use the first sounding with coverage above 10%:
	id = findfirst((df.VegCov .> 20.).&(df.VegCov .< 30.))
	
	Plots.plot(vegPolygon, legend=false, fill = (34,0.5,:gray))
	Sdng = getSdngVerts(df,id)
	testPoints = getPoints(Sdng , UInt8(10))
	#Plots.scatter!(testPoints)
	
	touchPoly = getPoly(Sdng.lon,Sdng.lat)
	Plots.plot!(touchPoly, linewidth=4)
	title!("Coverage: "*string(df.VegCov[id])*"%")
end

# ╔═╡ 620677cc-6817-11eb-0205-e9daaf8aaa37
begin
let
	## use just one month of data for faster gridding:
	idxTime = findall(month.(df.Date) .== 5)
	## or two month?:
	#in_month = in([4,5,6])
	#idxTime = findall(in_month.(month.(df.Date)).==true)
	## Or half a month?:
	#idxTime = findall((month.(df.Date) .== 5) .& (day.(df.Date) .<= 15))

	out, lons, lats = oversample(df[idxTime,:],tryparse(Float64,resolution), 										 mbr,"sif")
	
		
	Plots.heatmap(lons, lats, out, color= :viridis)#, clim=(-0.2,0.9))
	#fin = oversample(df, resolution, mbr,"sif")
	#Plots.plot(fin, color= :viridis)
	Plots.plot!(vegPolygon, legend=false, linewidth=3)
	Plots.plot!(refPolygon, legend=false, linewidth=3, linecolor=:blue)
	title!("n="*string(length(idxTime)))
end
end

# ╔═╡ f1f9738c-669b-11eb-3a6f-011326f23714
md"## Uncertainty/measurement noise:"

# ╔═╡ d1cec2be-71a7-11eb-1d36-e9c097be9b7f
md"
> **Negative SIF values are plausible due to the measurement uncertainty. Filtering negative SIF values will result in a positive bias!**
"

# ╔═╡ a9471dee-6b0a-11eb-0dde-19829d25aab5
begin

	#rnorm(n, mu, sig) = rand(Normal(mu, sig), n)
	
	### what can be expected (mean error estimate, uncertainty prediction):
	meanErr = mean(df[indRef,"sif_err"])
	
	### what is actually observed:
	fitErr = fit(Normal, df[indRef,"sif"])

	## plot:
	histogram(df[indRef,"sif"], bins = 100, normalize = :pdf, label="Observed σ = "*string(round(fitErr.σ; digits=2))*", μ = "*string(round(fitErr.μ; digits=2)))
	plot!(Normal(0.,meanErr), linewidth=4, label="Predicted σ = "*string(round(meanErr; digits=2))*", μ = 0")	
	title!("n="*string(length(indRef)))

end

# ╔═╡ Cell order:
# ╟─6659c6d4-7182-11eb-3b65-d939a40b4325
# ╟─1c2705ee-61a4-11eb-33f9-835cc5caa8a8
# ╟─b7921386-7186-11eb-20d2-d973f7940a35
# ╟─e7799730-7187-11eb-02a8-e7f76ed6081c
# ╟─b434a196-61a4-11eb-19fe-976e23dec488
# ╠═c25b5f6c-61a4-11eb-2c34-cbf860b5d314
# ╟─a4854224-7184-11eb-291a-e3667b3c33e2
# ╠═5250d5d2-718f-11eb-03ab-5d29ef4ce6ad
# ╟─34be05d8-7185-11eb-3bd3-d1e24885ed75
# ╟─8cc541f6-7185-11eb-3034-675a8bf6ff7c
# ╟─2be873b4-718d-11eb-1a33-a742dd4610ff
# ╟─bb5cb8a2-718d-11eb-3694-cbe3adf78a8d
# ╟─a33e9ac2-718c-11eb-18d4-b90c627f52ec
# ╟─8f48bb84-718e-11eb-333e-8f1f2df863d9
# ╟─38f5b8de-7190-11eb-0437-497aa85610f0
# ╟─2df3dee4-61a4-11eb-2003-5b56cfb49d01
# ╠═5b5b38b4-61a4-11eb-284b-e329267e0fb2
# ╟─4f00b9e0-61a4-11eb-0d33-f3ec82cd5396
# ╠═2f0d9a3a-61a5-11eb-2fc1-650e465800c8
# ╟─62fa453c-61a5-11eb-0b39-015a9ff4155c
# ╠═a18e5c16-61c3-11eb-37cf-edd1266f9dde
# ╠═179972c6-7257-11eb-3a3d-45c672612fb3
# ╟─bc20f27e-61d1-11eb-2497-179d53ab2360
# ╠═dfaee354-61d1-11eb-2a88-8bb2002c2136
# ╟─a8dc819e-7188-11eb-2c1b-8998b48d99f5
# ╟─cfd6f918-7189-11eb-33f4-c5133381580b
# ╟─c475fc80-718a-11eb-2c79-3f45c09aaca4
# ╟─2680b0a4-718c-11eb-2bad-137068539830
# ╟─471fe14a-718c-11eb-23a2-a3c5ed97e288
# ╟─a748c364-71a2-11eb-21fb-2ddc00692957
# ╟─a8564f6a-71a2-11eb-1797-67abf06ac8c1
# ╟─d753c838-71a2-11eb-151d-6394549ab869
# ╟─d7f995ce-71a2-11eb-3e35-b113a803a9fb
# ╟─91b73620-7193-11eb-070f-573411f9e611
# ╠═df1138ee-61d5-11eb-34b6-67ce049efb74
# ╠═14e5f60c-6687-11eb-37f6-fbf0d9061530
# ╠═ad814292-65d9-11eb-1390-35f78905835c
# ╠═e083f3fc-65da-11eb-290c-8df8464a4ff1
# ╟─53d0a158-62a6-11eb-3f8b-af81223169ec
# ╟─213784d6-7198-11eb-16a2-254812fbed1d
# ╟─b88cdc64-6509-11eb-24c7-fd78be2de33c
# ╠═80cc3a5c-629f-11eb-3807-f51439287448
# ╠═b1b2bb48-64fd-11eb-1ce8-274175a2d8ff
# ╟─b3fa469c-61d4-11eb-17ca-31a9319d9560
# ╠═50938100-62af-11eb-0027-2bc902261db2
# ╟─6288c4c6-6585-11eb-24c3-4b39125563bb
# ╠═35ce9bee-62a7-11eb-0118-9528ab9f33d8
# ╠═254018d0-62b7-11eb-2b1d-138ad8d088b1
# ╠═82ade7fe-62b7-11eb-3aa4-e9b8ba1f5a58
# ╟─50f46e82-64c3-11eb-29d9-933218ab2e66
# ╟─b30830b2-71a3-11eb-2699-53e1f7bcdb43
# ╟─a33e67bc-650a-11eb-2ad2-43fe662d8a78
# ╠═00f8867e-64f3-11eb-18f0-2f31200f7628
# ╠═7ff0299c-7265-11eb-326b-e5659d2f2f58
# ╠═7563ebee-64f5-11eb-30ca-f5e9d7856dd7
# ╟─eea1155e-64f5-11eb-3eed-f34a6000d5ea
# ╟─a3b5d552-650e-11eb-0db6-cdd12b877af8
# ╟─c21d0ec8-71a4-11eb-1aa7-b1f03c83a3b7
# ╠═a12da5e2-71a4-11eb-3f91-437744720d0d
# ╠═eb55f0c2-6583-11eb-2572-1b96530ce77a
# ╟─f62cce7c-725a-11eb-3b3b-ad594257cda6
# ╠═f33daa04-725a-11eb-2831-6fc756438afe
# ╠═724200f0-7262-11eb-1b6a-3723b249389d
# ╠═f3550094-71a5-11eb-1bd6-8150f1ad5604
# ╠═b99cff74-650e-11eb-3373-1d52b025e9a5
# ╟─4a8de51e-65db-11eb-2e99-c70cb7193464
# ╟─26556eb4-664f-11eb-16f6-2b8c0c1a5945
# ╠═4a72d39c-65da-11eb-07c2-2dd1c2bcac1e
# ╟─ac7223a8-668d-11eb-1906-a51e5b2ce8c1
# ╠═9555b186-664d-11eb-088a-bd3beb009eab
# ╟─c88c215e-6650-11eb-35fe-4dbe87cfb463
# ╟─a5b67f32-668f-11eb-090f-07b24fe0660d
# ╠═964d3d38-665d-11eb-1eb4-adbb37c289b6
# ╟─98684c54-6660-11eb-175f-c96664cc2f04
# ╟─0ff6e15e-6661-11eb-3e48-3b73b977b0fb
# ╟─0cda4bd4-6691-11eb-20c7-9b5097e1e7e9
# ╟─388e8398-6692-11eb-16c8-2bcccd929b61
# ╠═9ba7d79a-6692-11eb-30e8-8d548156f71c
# ╟─46c0b6a0-6824-11eb-2407-a37dca8ef4bf
# ╠═6988f26c-682c-11eb-3261-ed8a7c839fb4
# ╠═2140d9a2-682c-11eb-0be0-8d4cba6681bb
# ╠═5623b812-6825-11eb-091d-3b676411207c
# ╠═d51c868a-760d-11eb-3420-9b51e5d9299a
# ╠═e28b9fc0-6a57-11eb-183f-6f47f1c39d6a
# ╟─e53f342c-669a-11eb-3afb-fbeede548218
# ╟─59c35358-682d-11eb-0fd2-77744637cdbe
# ╠═620677cc-6817-11eb-0205-e9daaf8aaa37
# ╠═2e2828e8-67ee-11eb-1b37-3ffc4a162aeb
# ╟─f1f9738c-669b-11eb-3a6f-011326f23714
# ╟─d1cec2be-71a7-11eb-1d36-e9c097be9b7f
# ╠═a9471dee-6b0a-11eb-0dde-19829d25aab5
