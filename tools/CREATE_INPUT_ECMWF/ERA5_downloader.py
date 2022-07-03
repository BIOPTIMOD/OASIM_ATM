import cdsapi

c = cdsapi.Client()

year='2019'
month=['01','02', '03','04', '05', '06','07', '08', '09','10', '11', '12']
month=['01']
day = [ '%02d' %f for f in range(1,2)] #32
time =[ '%02d:00' %f for f in range(24)]

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': ['total_column_cloud_liquid_water', 'total_column_ozone'],
        'year': year,
        'month': month,
        'day': day,
        'time': time,
        'area': [47, -20, 29, 40],
    },
    'ERA5_Med.nc')