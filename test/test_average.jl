using XRadio: Testing, subdataset

@testset "average" begin
    t0 = 1.6148160e9
    # Two scans of two 10 s samples each, the second scan starting 100 s in.
    function two_scans(; eltype = ComplexF32)
        m = Testing.measurement_set(;
            times = t0 .+ [0.0, 10.0, 100.0, 110.0], integration_time = 10.0, eltype
        )
        m = copy(m)
        m[:scan_name] = rebuild(m[:scan_name], ["1", "1", "2", "2"])
        V = copy(m[:visibility])
        W = copy(m[:weight])
        F = copy(m[:flag])
        U = copy(m[:uvw])
        cell(f, b, t) = (Polarization(At("RR")), Frequency(f), BaselineID(At(b)), Ti(t))
        # Scan 1, first channel: a weighted mean, one sample flagged, both flagged.
        V[cell(1, 0, 1)...] = 1 + 1im
        V[cell(1, 0, 2)...] = 3 + 3im
        W[cell(1, 0, 2)...] = 3
        F[cell(1, 1, 1)...] = true
        V[cell(1, 1, 2)...] = 5
        W[cell(1, 1, 2)...] = 2
        F[cell(1, 2, 1)...] = true
        F[cell(1, 2, 2)...] = true
        # Second channel: a sample of zero weight, a sample that is not finite.
        V[cell(2, 0, 1)...] = 100
        W[cell(2, 0, 1)...] = 0
        V[cell(2, 0, 2)...] = 2
        V[cell(2, 1, 1)...] = NaN
        V[cell(2, 1, 2)...] = 4
        # Scan 2: every datum of its second sample weighs three times its first.
        view(W, Ti(4)) .= 3
        U[BaselineID(At(0)), Ti(1)] = [1, 2, 3]
        U[BaselineID(At(0)), Ti(2)] = [4, 5, 6]
        U[BaselineID(At(0)), Ti(4)] = [4, 4, 4]
        m[:visibility] = V
        m[:weight] = W
        m[:flag] = F
        m[:uvw] = U
        return m
    end
    ms = two_scans()
    out = average(ms, ByScan())
    at(A, f, b, t) = A[Polarization(At("RR")), Frequency(f), BaselineID(At(b)), Ti(t)]

    @testset "correlated data" begin
        V, W, F = out[:visibility], out[:weight], out[:flag]
        @test size(V) == (4, 8, 3, 2)
        @test eltype(V) == ComplexF32
        @test at(V, 1, 0, 1) ≈ 2.5 + 2.5im
        @test at(W, 1, 0, 1) == 4
        @test !at(F, 1, 0, 1)
        @test (at(V, 1, 1, 1), at(W, 1, 1, 1), at(F, 1, 1, 1)) == (5, 2, false)
        @test isnan(real(at(V, 1, 2, 1))) && isnan(imag(at(V, 1, 2, 1)))
        @test (at(W, 1, 2, 1), at(F, 1, 2, 1)) == (0, true)
        @test (at(V, 2, 0, 1), at(W, 2, 0, 1)) == (2, 1)
        @test (at(V, 2, 1, 1), at(W, 2, 1, 1)) == (4, 1)
        @test at(V, 3, 1, 1) == 1 && at(W, 3, 1, 1) == 2
        @test all(==(4), W[Ti(2)]) && all(==(1), V[Ti(2)]) && !any(F[Ti(2)])
    end

    @testset "uvw, epochs and integration time" begin
        # Weighted per baseline and time by the weight contributing over every
        # polarization and channel: 31 at the first sample, 34 at the second.
        @test out[:uvw][BaselineID(At(0)), Ti(1)] ≈ (31 .* [1, 2, 3] .+ 34 .* [4, 5, 6]) ./ 65
        @test out[:uvw][BaselineID(At(0)), Ti(2)] ≈ [3, 3, 3]
        @test all(iszero, out[:uvw][BaselineID(At(1))])
        # The mid-points of the spans, whatever the weights.
        @test times(out) ≈ [t0 + 5, t0 + 105]
        it = metadata(lookup(out, Ti))[:integration_time]
        @test it isa Measure && XRadio.value(it) == 20.0
        @test metadata(it)[:units] == "s"
        @test metadata(lookup(out, Ti))[:scale] == "utc"
    end

    @testset "per-time labels" begin
        @test collect(out[:scan_name]) == ["1", "2"]
        @test collect(out[:field_name]) == ["field_0", "field_0"]
        @test_throws "`scan_name` takes more than one value under the label :all" average(
            ms, _ -> :all
        )
    end

    @testset "one span per Measurement Set" begin
        @test_throws "one Measurement Set states one `integration_time`" average(
            ms, _ -> [1, 1, 1, 2]
        )
    end

    @testset "interleaved labels" begin
        m = Testing.measurement_set(; times = t0 .+ (0:10:30))
        @test times(average(m, _ -> [1, 2, 1, 2])) ≈ [t0 + 10, t0 + 20]
    end

    @testset "epochs do not depend on weights or flags" begin
        m = copy(Testing.measurement_set(; times = t0 .+ (0:10:30)))
        W = copy(m[:weight])
        view(W, Ti(1)) .= 0
        view(W, Ti(4)) .= 5
        F = copy(m[:flag])
        view(F, Ti(2)) .= true
        m[:weight] = W
        m[:flag] = F
        plain = Testing.measurement_set(; times = t0 .+ (0:10:30))
        @test times(average(m, _ -> :all)) == times(average(plain, _ -> :all))
    end

    @testset "a label for the whole Measurement Set" begin
        one = average(Testing.measurement_set(; times = t0 .+ (0:10:30)), _ -> :all)
        @test length(times(one)) == 1
        @test times(one) ≈ [t0 + 15]
        @test XRadio.value(metadata(lookup(one, Ti))[:integration_time]) == 40
    end

    @testset "conforms, and keeps what is not stated per time" begin
        @test isempty(check(out))
        @test eltype(average(two_scans(; eltype = ComplexF64), ByScan())[:visibility]) == ComplexF64
        @test out[:baseline_antenna1_name] == ms[:baseline_antenna1_name]
        @test metadata(out)[:data_groups] == metadata(ms)[:data_groups]
        @test XRadio._is_linked(branches(out)[:field_and_source_base], out)
        @test ByScan()(out) == unique(ByScan()(ms))
        @test parent(ms[:visibility])[1, 1, 1, 1] == 1 + 1im
    end

    @testset "sub-datasets on clocks of their own" begin
        weather = subdataset(
            :weather;
            station_name = ["A1"], cartesian_pos_label = ["x", "y", "z"],
            time_weather = [t0, t0 + 300], station_position = reshape([1.0, 2.0, 3.0], 3, 1),
            temperature = reshape([283.15, 284.15], 2, 1),
        )
        withweather = copy(ms)
        withweather[:weather] = weather
        got = average(withweather, ByScan())
        @test isempty(check(got))
        @test parent(got.weather[:temperature]) == parent(weather[:temperature])
        @test lookup(got.weather, TimeWeather) == lookup(weather, TimeWeather)
        @test XRadio._is_linked(branches(got)[:weather], got)
    end

    @testset "a linked sub-dataset stated per time" begin
        linked = copy(ms)
        extra = XRadio.Dataset(; metadata = Dict{Symbol, Any}(:type => "extra"))
        extra[:observed] = DimArray(["a", "a", "b", "b"], (dims(ms, Ti),))
        linked[:extra] = extra
        got = average(linked, ByScan())
        @test XRadio._is_linked(branches(got)[:extra], got)
        @test collect(got.extra[:observed]) == ["a", "b"]
        extra[:direction] = DimArray(
            [1.0 1.0 0.0 4.0; 2.0 2.0 0.0 4.0],
            (SkyDirLabel(["ra", "dec"]), dims(ms, Ti))
        )
        linked[:extra] = extra
        @test_throws "`extra/direction` is stated per time" average(linked, ByScan())
    end

    @testset "variables it does not average" begin
        bt = (dims(ms, BaselineID), dims(ms, Ti))
        full = copy(ms)
        full[:effective_integration_time] = DimArray(
            fill(10.0, 3, 4), bt; metadata = Dict{Symbol, Any}(:units => "s", :type => "quantity")
        )
        @test_throws "`effective_integration_time` is stated per time" average(full, ByScan())
        g = copy(ms)
        g[:visibility_corrected] = rebuild(ms[:visibility], 2 .* parent(ms[:visibility]))
        add_datagroup!(g, :corrected; correlated_data = "VISIBILITY_CORRECTED")
        @test_throws "one data group" average(g, ByScan())
    end

    @testset "an unknown variable stated per time" begin
        odd = copy(ms)
        odd[:mystery] = DimArray(rand(4), (dims(ms, Ti),))
        @test_throws "`mystery` is stated per time, and `average` does not know how to average it" average(
            odd, ByScan()
        )
    end

    @testset "a processing set" begin
        ps = Testing.processing_set(; nscan = 2, nspw = 2, ntime = 3)
        avg = average(ps, ByScan())
        @test collect(keys(avg)) == collect(keys(ps))
        @test collect(keys(groupby(avg, ByScan()))) == collect(keys(groupby(ps, ByScan())))
        @test all(ms -> length(times(ms)) == 1, avg)
        @test isempty(check(avg))
    end
    @testset "over channels" begin
        f = frequencies(ms)
        width = XRadio.value(metadata(lookup(ms, Frequency))[:channel_width])
        two = average(ms, ChannelBins(2))
        V, W, F = two[:visibility], two[:weight], two[:flag]
        @test size(V) == (4, 4, 3, 4)
        @test frequencies(two) == [(f[i] + f[i + 1]) / 2 for i in 1:2:8]
        @test XRadio.value(metadata(lookup(two, Frequency))[:channel_width]) == 2 * width
        @test times(two) == times(ms)
        @test parent(two[:uvw]) == parent(ms[:uvw])
        @test (at(V, 1, 0, 1), at(W, 1, 0, 1), at(F, 1, 0, 1)) == (1 + 1im, 1, false)
        @test isnan(real(at(V, 1, 1, 1))) && (at(W, 1, 1, 1), at(F, 1, 1, 1)) == (0, true)
        @test at(V, 1, 0, 2) ≈ (11 + 9im) / 4 && at(W, 1, 0, 2) == 4
        @test isempty(check(two))
    end

    @testset "time then channels equals both at once" begin
        both = average(ms, ByScan(), ChannelBins(2))
        @test isequal(parent(both[:visibility]), parent(average(average(ms, ByScan()), ChannelBins(2))[:visibility]))
        @test isequal(parent(average(ms, ChannelBins(2), ByScan())[:visibility]), parent(both[:visibility]))
        v, w, fl = parent(ms[:visibility]), parent(ms[:weight]), parent(ms[:flag])
        scan_of = [1, 1, 2, 2]
        for I in CartesianIndices(parent(both[:visibility]))
            p, c, b, s = Tuple(I)
            num, den = zero(ComplexF64), 0.0
            for ch in (2c - 1):(2c), t in findall(==(s), scan_of)
                wi, vi = w[p, ch, b, t], v[p, ch, b, t]
                (fl[p, ch, b, t] || !(wi > 0) || !isfinite(vi)) && continue
                num += wi * vi
                den += wi
            end
            @test parent(both[:weight])[I] == den
            iszero(den) || @test parent(both[:visibility])[I] ≈ num / den
        end
        @test isempty(check(both))
    end

    @testset "channel labels" begin
        apart = average(ms, ByChannel(_ -> [1, 2, 1, 2, 3, 4, 3, 4]))
        f = frequencies(ms)
        @test frequencies(apart) == [(f[1] + f[3]) / 2, (f[2] + f[4]) / 2, (f[5] + f[7]) / 2, (f[6] + f[8]) / 2]
        @test length(frequencies(average(ms, ByChannel(_ -> :all)))) == 1
        @test XRadio.averaging_axis(ByScan()) == :time
        @test XRadio.averaging_axis(ChannelBins(2)) == :frequency
        @test_throws "positive number of channels" ChannelBins(0)
        @test_throws "gathers 3 channels" average(ms, ChannelBins(3))
        @test_throws "same frequency" average(ms, ByChannel(_ -> [1, 2, 2, 1, 3, 3, 4, 4]))
        @test_throws "one label per channel" average(ms, ByChannel(_ -> [1, 2]))
        odd = copy(ms)
        odd[:mystery] = DimArray(rand(8), (dims(ms, Frequency),))
        @test_throws "`mystery` is stated per frequency" average(odd, ChannelBins(2))
        ps = Testing.processing_set(; nscan = 2, nspw = 2, ntime = 3)
        @test isempty(check(average(ps, ByScan(), ChannelBins(4))))
    end
end
