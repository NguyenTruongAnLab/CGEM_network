"""Analyze validation data from March 2025 field campaign"""
import pandas as pd
import numpy as np

# Load validation data
df = pd.read_csv('INPUT/Cases/Mekong_Delta_Full/Validation/Mekong_4_branches_Mar_2025.csv')

print('=== VALIDATION DATA SUMMARY (March 2025) ===\n')

for river in ['Hau_River', 'Ham_Luong', 'Co_Chien', 'My_Tho']:
    rdf = df[df['River'] == river].sort_values('Distance')
    dist_min = rdf['Distance'].min()
    dist_max = rdf['Distance'].max()
    print(f'\n--- {river} ({len(rdf)} stations, {dist_min:.1f}-{dist_max:.1f} km) ---')
    
    # Salinity
    sal = rdf['Surface Salinity'].dropna()
    print(f'  Salinity: {sal.min():.1f} - {sal.max():.1f} PSU')
    
    # O2 (mg/L -> umol/L: multiply by 31.25)
    o2 = rdf['O2 (mg/L)'].dropna()
    print(f'  O2: {o2.min():.2f} - {o2.max():.2f} mg/L ({o2.min()*31.25:.0f} - {o2.max()*31.25:.0f} umol/L)')
    
    # TSS
    tss = rdf['TSS (mg/L)'].dropna()
    if len(tss) > 0:
        print(f'  TSS: {tss.min():.1f} - {tss.max():.1f} mg/L')
    
    # pCO2
    pco2 = rdf['pCO2 (ppm)'].dropna()
    print(f'  pCO2: {pco2.min():.0f} - {pco2.max():.0f} ppm')
    
    # pH
    ph = rdf['pH'].dropna()
    print(f'  pH: {ph.min():.2f} - {ph.max():.2f}')
    
    # Alkalinity
    alk = rdf['Alkalinity (mmol/L)'].dropna()
    if len(alk) > 0:
        print(f'  Alkalinity: {alk.min():.2f} - {alk.max():.2f} mmol/L ({alk.min()*1000:.0f} - {alk.max()*1000:.0f} umol/L)')
    
    # Nutrients (mgN/L -> umol/L: multiply by 71.4)
    no3 = rdf['NO3-N (mgN/L)'].dropna()
    if len(no3) > 0:
        print(f'  NO3: {no3.min():.2f} - {no3.max():.2f} mgN/L ({no3.min()*71.4:.1f} - {no3.max()*71.4:.1f} umolN/L)')
    
    nh4 = rdf['NH4 (mgN/L)'].dropna()
    if len(nh4) > 0:
        print(f'  NH4: {nh4.min():.2f} - {nh4.max():.2f} mgN/L ({nh4.min()*71.4:.1f} - {nh4.max()*71.4:.1f} umolN/L)')
    
    # DOC (mgC/L -> umol/L: multiply by 83.3)
    doc = rdf['DOC (mgC/L)'].dropna()
    if len(doc) > 0:
        print(f'  DOC: {doc.min():.2f} - {doc.max():.2f} mgC/L ({doc.min()*83.3:.0f} - {doc.max()*83.3:.0f} umolC/L)')
    
    # Chl-a
    chla = rdf['Chl-a (ug/L)'].dropna()
    if len(chla) > 0:
        print(f'  Chl-a: {chla.min():.2f} - {chla.max():.2f} ug/L')
    
    # CH4 (ugC/L -> nmol/L: divide by 12 * 1000)
    ch4 = rdf['ugC-CH4/L'].dropna()
    if len(ch4) > 0:
        ch4_nmol = ch4 / 12 * 1000  # ugC -> nmol CH4
        print(f'  CH4: {ch4.min():.2f} - {ch4.max():.2f} ugC/L ({ch4_nmol.min():.0f} - {ch4_nmol.max():.0f} nmol/L)')
    
    # N2O (ugN/L -> nmol/L: divide by 28 * 1000)
    n2o = rdf['ugN-N2O/L'].dropna()
    if len(n2o) > 0:
        n2o_nmol = n2o / 28 * 1000  # ugN -> nmol N2O
        print(f'  N2O: {n2o.min():.2f} - {n2o.max():.2f} ugN/L ({n2o_nmol.min():.0f} - {n2o_nmol.max():.0f} nmol/L)')

print('\n\n=== KEY OBSERVATIONS FOR MODEL CALIBRATION ===')
print('''
1. SALINITY INTRUSION:
   - Ocean salinity ~30 PSU at mouth
   - Hau_River: drops to 0.1 PSU by ~58 km (intrusion ~40-50 km)
   - Ham_Luong: drops to 0.1 PSU by ~66 km (intrusion ~55-60 km)  
   - Co_Chien: drops to 0.1 PSU by ~70 km (intrusion ~65-70 km)
   - My_Tho: drops to 0.1 PSU by ~51 km (intrusion ~45-50 km)

2. OXYGEN:
   - HIGH at ocean (8+ mg/L = 250+ umol/L, ~100% saturation)
   - DROPS significantly upstream (5-6 mg/L = 160-190 umol/L, ~70-80%)
   - Minimum values in freshwater zone (Can Tho area)

3. pCO2:
   - LOW at ocean (~500 ppm, near equilibrium)
   - INCREASES dramatically upstream (3000-4700 ppm!)
   - Maximum near freshwater zone (organic decomposition)

4. TSS (SPM):
   - Variable 6-38 mg/L
   - Often HIGHER in mid-estuary (ETM zone) and freshwater
   - Lower at ocean mouth

5. NUTRIENTS (very low!):
   - NO3: 0.09-0.93 mgN/L (6-66 umol/L) - MUCH LOWER than model
   - NH4: 0.01-0.07 mgN/L (0.7-5 umol/L) - MUCH LOWER than model
   - PO4: 0.01-0.05 mgP/L (0.3-1.6 umol/L) - reasonable

6. DOC:
   - 1.3-2.7 mgC/L (108-225 umol/L) - LOWER than typical river

7. Chl-a:
   - 0.7-7 ug/L - reasonable for turbid estuary

8. CH4 and N2O:
   - CH4: 0.5-2.9 ugC/L (40-240 nmol/L)
   - N2O: 0.18-5.3 ugN/L (6-190 nmol/L)
''')
