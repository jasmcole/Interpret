function outregion = GetSavedRegion(region)

switch region
    case 'Synthetic FFT'
        outregion(1) = 341;
        outregion(2) = 137;
        outregion(3) = 21;
        outregion(4) = 27;
    case 'Synthetic'
        outregion(1) = 1;
        outregion(2) = 1;
        outregion(3) = 747;
        outregion(4) = 299;
    case 'Gemini 2012 20121114 Late Image'
        outregion(1) = 2069;
        outregion(2) = 1779;
        outregion(3) = 344;
        outregion(4) = 341;
    case 'Gemini 2012 20121114 Late FFT'
        outregion(1) = 160;
        outregion(2) = 154;
        outregion(3) = 8;
        outregion(4) = 35;
    case 'Gemini 2012 20121114 Image'
        outregion(1) = 454;
        outregion(2) = 1778;
        outregion(3) = 2178;
        outregion(4) = 313;
    case 'Gemini 2012 20121114 FFT'
        outregion(1) = 1014;
        outregion(2) = 141;
        outregion(3) = 52;
        outregion(4) = 32;
    case 'Gemini 2012 20121026 Early Image'
        outregion(1) = 2057;
        outregion(2) = 1707;
        outregion(3) = 377;
        outregion(4) = 213;
    case 'Gemini 2012 20121026 Early FFT'
        outregion(1) = 153;
        outregion(2) = 87;
        outregion(3) = 33;
        outregion(4) = 37;
    case 'Gemini 2012 20121026 Mid Image'
        outregion(1) = 2388;
        outregion(2) = 1700;
        outregion(3) = 447;
        outregion(4) = 227;
    case 'Gemini 2012 20121026 Mid FFT'
        outregion(1) = 192;
        outregion(2) = 95;
        outregion(3) = 27;
        outregion(4) = 43;
    case 'Gemini 2012 20121026 Late Image'
        outregion(1) = 2984;
        outregion(2) = 1675;
        outregion(3) = 377;
        outregion(4) = 230;
    case 'Gemini 2012 20121026 Late FFT'
        outregion(1) = 168;
        outregion(2) = 107;
        outregion(3) = 15;
        outregion(4) = 18;
    case 'none'
        outregion(1) = 0;
        outregion(2) = 0;
        outregion(3) = 0;
        outregion(4) = 0;
end

end