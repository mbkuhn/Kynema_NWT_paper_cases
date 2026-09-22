using DelimitedFiles
using Interpolations
using DataFrames
using LaTeXStrings
using StatsBase
using MAT
using Measurements
include("helper_scripts/FFTAnalysis.jl")
include("helper_scripts/get_evenly_spaced_time_series.jl")

# Define smoothing type and relevant parameters
# Options are: 
# Ensemble Averaging -> smoothing_type = 1; smoothing_parameter = 10000
# Constant Filter (MARIN's Method) -> smoothing_type = 2; smoothing_parameter = 30
# Gaussian -> smoothing_type = 3; smoothing_parameter = 30
smoothing_type = 2
smoothing_parameter = -1
if smoothing_type==1
    smoothing_parameter = 10000
else
    smoothing_parameter = 30
end

file_amrwind = "data/amr_wind_hos_nwt_sampling_time_series.txt"
file_amrwind2 = "data/amr_wind_hos_nwt_sampling_time_series_finer.txt"
file_hosnwt = "data/HOS_NWT_probes.dat"

rows_to_skip = 49
# Data from probes for JONSWAP spectrum from HOS-NWT simulation
HOS_NWT_Case4_data = readdlm(file_hosnwt)
HOS_NWT_Case4_data = HOS_NWT_Case4_data[rows_to_skip:end,1:3]
HOS_NWT_Case4_data = Float64.(HOS_NWT_Case4_data)
HOS_NWT_Case4_time = HOS_NWT_Case4_data[:,1]
HOS_NWT_Case4_wave = HOS_NWT_Case4_data[:,3]

# Rewrite HOS_NWT time due to lack of precision
HOS_NWT_dt = HOS_NWT_Case4_time[2] - HOS_NWT_Case4_time[1]
for ii in 1:length(HOS_NWT_Case4_time)
    HOS_NWT_Case4_time[ii] = HOS_NWT_Case4_time[1] + (ii-1)*HOS_NWT_dt
end

# AMR-Wind
amrwind_data = readdlm(file_amrwind)
amrwind_data = amrwind_data[2:end,1:4]
amrwind_data = Float64.(amrwind_data)
amrwind_time = amrwind_data[:,1]
amrwind_wave = amrwind_data[:,4]

amrwind_data2 = readdlm(file_amrwind2)
amrwind_data2 = amrwind_data2[2:end,1:4]
amrwind_data2 = Float64.(amrwind_data2)
amrwind_time2 = amrwind_data2[:,1]
amrwind_wave2 = amrwind_data2[:,4]

# Experiment
file_exp = "/Users/mkuhn/Library/CloudStorage/OneDrive-NLR/testruns_data/OC6_PhaseIa_dataset_and_more/OC6_LC33_LC53_CFD_scripts_data/LC33_FixedStructure/ZXP0_LC33_full.txt"
ZXP0_data = readdlm(file_exp)
ZXP0_data = ZXP0_data[2001:end,:]
ZXP0_data = Float64.(ZXP0_data)
ZXP0_data[:,1] = ZXP0_data[:,1].-200 # "." dot is to subtract elementwise for the column
EXP_time = ZXP0_data[:,1]
EXP_wave = ZXP0_data[:,2]

# Plot Waves Time series (Model Scale)
rescaling = 50.
Time_series_plot = plot(EXP_time./sqrt(rescaling),EXP_wave./(rescaling),label="Experiment",color=:black,legend=:topright)
plot!(HOS_NWT_Case4_time,HOS_NWT_Case4_wave,label="HOS-NWT", color=:blue)
plot!(amrwind_time,amrwind_wave,label="Kynema-SGF, \$\\Delta{z}\$ = 0.0188m", color=:green, linestyle=:dashdot)
plot!(amrwind_time2,amrwind_wave2,label="Kynema-SGF, \$\\Delta{z}\$ = 0.0094m", color=:red, linestyle=:dot)
xlabel!("Time [sec]")
ylabel!("Wave Elevation [m]")
savefig(Time_series_plot,"plotting_outputs/wave_elevation_time_series.pdf")

### **** PSD COMPUTATIONS AND PLOTS ******
stepRange = collect(3000:108000)
f, P1, PSD, fdouble, P2 = FFTAnalysis(EXP_time[stepRange]./sqrt(rescaling),EXP_wave[stepRange]./rescaling,true,false,smoothing_type,smoothing_parameter)
wave_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,ylims=(1e-7,1e-2),label="Experiment",xlims=(0,1.5),color=:black,linewidth=2,legend=:bottomright)

f, P1, PSD, fdouble, P2 = FFTAnalysis(HOS_NWT_Case4_time,HOS_NWT_Case4_wave,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="HOS-NWT",color=:blue,linewidth=2)

# Coarser AMR-Wind simulation
dt_array = zeros(length(amrwind_data[:,1]))
for ii=1:length(dt_array)-1
    local dt = amrwind_data[ii+1,1].-amrwind_data[ii,1]
    dt_array[ii] = dt
end
t_spacing = median(dt_array)

AMRWind_WAVE_time_interp = Interpolations.deduplicate_knots!(amrwind_time)
itp_waves_amrwind =  interpolate((AMRWind_WAVE_time_interp,), amrwind_wave, Gridded(Linear()))
amr_time_evenly_spaced = collect(amrwind_data[1,1]:t_spacing:amrwind_data[end,1])
amrwind_wave_evenly_spaced = itp_waves_amrwind.(amr_time_evenly_spaced)
AMRWind_time = amr_time_evenly_spaced
AMRWind_Wave = amrwind_wave_evenly_spaced

f, P1, PSD, fdouble, P2 = FFTAnalysis(AMRWind_time,AMRWind_Wave,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Kynema-SGF, \$\\Delta{z}\$ = 0.0188m",color=:green,linestyle=:dashdot,linewidth=2)

# Finer AMR-Wind simulation
dt_array = zeros(length(amrwind_data2[:,1]))
for ii=1:length(dt_array)-1
    local dt = amrwind_data2[ii+1,1].-amrwind_data2[ii,1]
    dt_array[ii] = dt
end
t_spacing = median(dt_array)

AMRWind_WAVE_time_interp = Interpolations.deduplicate_knots!(amrwind_time2)
itp_waves_amrwind =  interpolate((AMRWind_WAVE_time_interp,), amrwind_wave2, Gridded(Linear()))
amr_time_evenly_spaced = collect(amrwind_data2[1,1]:t_spacing:amrwind_data2[end,1])
amrwind_wave_evenly_spaced = itp_waves_amrwind.(amr_time_evenly_spaced)
AMRWind_time = amr_time_evenly_spaced
AMRWind_Wave = amrwind_wave_evenly_spaced

f, P1, PSD, fdouble, P2 = FFTAnalysis(AMRWind_time,AMRWind_Wave,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Kynema-SGF, \$\\Delta{z}\$ = 0.0094m",color=:red,linestyle=:dot,linewidth=2)

xlabel!("Frequency [Hz]")
ylabel!("PSD of Incident Wave [m\$^2\$/Hz]")
savefig(wave_plot,"plotting_outputs/wave_PSD.pdf")