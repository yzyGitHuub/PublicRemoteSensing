// Copyright © 2020 Zhaoyuan Yao, All right reserved.
//
var ndvi_min = 0.1;
var ndvi_max = 0.9;

// daman
var para = ee.Feature(ee.Geometry.Point([100.37223,38.85551]), {'top_k': 17.569135851448095, 'top_b': 2.2476314084363076, 'bottom_k': 3.7577242998718323, 'bottom_b': 0.5151560689543699})

var wet_k = ee.Image(ee.Number(para.get('top_k')));
var wet_b = ee.Image(ee.Number(para.get('top_b')));
var dry_k = ee.Image(ee.Number(para.get('bottom_k')));
var dry_b = ee.Image(ee.Number(para.get('bottom_b')));

var filterDate = ee.Filter.date('2017-06-05', '2020-09-01');
// print(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(filterDate))
var image = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
              .filter(ee.Filter.bounds(para.geometry()))
              .filter(ee.Filter.lte('CLOUD_COVER', 3))
              .sort('system:time_start')
              .filter(filterDate).first();
var dilatedCloudBitMask = (1 << 1);
var cirrusBitMask = (1 << 2);
var cloudBitMask = (1 << 3);
var cloudShadowBitMask = (1 << 4);
var snowBitMask = (1 << 5);
var waterBitMask = (1 << 7);
var qa = image.select('QA_PIXEL');
var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(snowBitMask).eq(0))
      .and(qa.bitwiseAnd(waterBitMask).eq(0))
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0))
      .and(qa.bitwiseAnd(dilatedCloudBitMask).eq(0))
      .and(qa.bitwiseAnd(cloudBitMask).eq(0));
var date = ee.Date(ee.Image(image).get('system:time_start'));
print(date, image)

var red = ee.Image(image).select('SR_B4').multiply(2.75e-05).add(-0.2);
var nir = ee.Image(image).select('SR_B5').multiply(2.75e-05).add(-0.2);
var swir = ee.Image(image).select(['SR_B7']).multiply(2.75e-05).add(-0.2);

var ndvi = nir.subtract(red).divide(nir.add(red));
// ndvi = ndvi.multiply(ndvi.gte(ndvi_min)).add(ndvi.lt(ndvi_min).multiply(ndvi_min))
// ndvi = ndvi.multiply(ndvi.lte(ndvi_max)).add(ndvi.gt(ndvi_max).multiply(ndvi_max))
var fr = ndvi.subtract(ndvi_min).divide(ndvi_max - ndvi_min).pow(1);
var str  = swir.subtract(1).multiply(swir.subtract(1)).divide(swir.multiply(2));

var phi_min = fr.multiply(1.26);
var optram = ndvi.multiply(dry_k).add(dry_b)
                  .subtract(str)
                  .divide(dry_b.subtract(wet_b)
                               .add(dry_k.subtract(wet_k).multiply(ndvi))// fr
                  )
optram = optram.multiply(optram.gt(0))
optram = optram.multiply(optram.lte(1)).add(optram.gt(1).multiply(1))
var phi = optram.multiply(ee.Image(1.26).subtract(phi_min))
                .add(phi_min);

var mo_image = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
                 .filter(ee.Filter.date(date, date.advance(1, 'day')))
                 .first();
var p = mo_image.select('surface_pressure').divide(1000) // Pa -> kPa
// ACTUAL VAPOR PRESSURE [KPA]
var tdp = mo_image.select('dewpoint_temperature_2m')
var ea = tdp.expression(
    '0.611 * pow(10, 7.45 / (235 + Td))',{
    'Td': tdp})
// HEAT CAPACITY OF MOIST AIT AT CONSTANT PRESSURE [J KG-1 K-1]
var q = ee.Image(ea).multiply(0.622).divide(ee.Image(p).add(ee.Image(ea).multiply(0.622 - 1)))
var cp = ee.Image(1).subtract(q).multiply(1003.5).add(ee.Image(q).multiply(1865))
// SLOPE OF SATURATED VAPOR PRESWSURE VERSUS AIR TEMPERATURE [KPA K-1]
var t_air_k = mo_image.select('temperature_2m')
var delta = t_air_k.expression(
  '4098.0 * (0.6108 * exp(17.27 * T_C / (T_C + 237.3))) / ((T_C + 237.3) ** 2)', {
    'T_C': t_air_k.subtract(273.15)
  }).resample('bilinear')
// LATENT HEAT OF VAPORISATION [J KG-1]
var lambda = ee.Image(t_air_k).subtract(273.15).multiply(-2.361e-3).add(2.501).multiply(1e6)
// PSICROMETRIC CONSTANT [KPA K-1]
var gamma = ee.Image(cp).multiply(p).divide(ee.Image(lambda).multiply(0.622)).resample('bilinear')
// NET RADIATION [W M-2]
var Rn = mo_image.select('surface_net_solar_radiation_sum').add(mo_image.select('surface_net_thermal_radiation_sum')).multiply(1.1574E-5).resample('bilinear')

var EF = ee.Image(phi).multiply(delta).divide(delta.add(gamma))

// // v: 0.05, s: 0.4
// var G  = fr.multiply(-1).add(1).multiply(0.4 - 0.05).add(0.05).multiply(Rn)
// var LE = Rn.subtract(G).multiply(EF); // latent heat flux[W m-2]
// var H  = Rn.subtract(G).multiply(EF.multiply(-1).add(1)); // sensible heat flux[W m-2]

var LE = Rn.multiply(EF).rename('LE'); // latent heat flux[W m-2]
var H  = Rn.multiply(EF.multiply(-1).add(1)); // sensible heat flux[W m-2]

var visualization = {
  min: 0.0,
  max: 5,
  palette: [
    'a50026', 'd73027', 'f46d43', 'fdae61', 'fee08b', 'ffffbf',
    'd9ef8b', 'a6d96a', '66bd63', '1a9850', '006837',
  ]
};
                            
Map.addLayer(LE.updateMask(mask).divide(28.94), visualization, 'OPTRAM-ET')
Map.addLayer(para)
Map.addLayer(ee.Image('MODIS/006/MOD16A2/2017_07_12').select('LE')
              .multiply(10000)
              .multiply(1.1574e-5)
              .divide(28.94)
              , visualization, 'MOD16')
Map.addLayer(ee.Image('CAS/IGSNRR/PML/V2_v017/2017-07-12').select('Ei')
        .add(ee.Image('CAS/IGSNRR/PML/V2_v017/2017-07-12').select('Ec'))
        .add(ee.Image('CAS/IGSNRR/PML/V2_v017/2017-07-12').select('Es'))
        , visualization, 'PML')

var export_region = ee.Geometry.BBox(100.5440851, 38.7009498, 100.6239076 , 38.7526427)
Map.addLayer(export_region)
var station_name = 'Daman_small';

var ET_AOPTRAMET = LE.updateMask(mask).divide(28.94)
Export.image.toDrive({
  scale : 30,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: ET_AOPTRAMET,
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_LC08_133033_20170716',
})

var ET_MOD16A2 = ee.Image('MODIS/006/MOD16A2/2017_07_12').select('LE')
                   .multiply(10000)
                   .multiply(1.1574e-5)
                   .divide(28.94)
Export.image.toDrive({
  scale : 500,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: ET_MOD16A2,
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_MOD16A2_2017_07_12',
})

var ET_PMLv2 = ee.Image('CAS/IGSNRR/PML/V2_v017/2017-07-12').select('Ei')
          .add(ee.Image('CAS/IGSNRR/PML/V2_v017/2017-07-12').select('Ec'))
          .add(ee.Image('CAS/IGSNRR/PML/V2_v017/2017-07-12').select('Es'))
Export.image.toDrive({
  scale : 500,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: ET_PMLv2,
  folder: 'global_et_example',
  description: station_name,
  fileNamePrefix: station_name + '_PMLv2_2017-07-12',
})

var TrueColor = image.select(['SR_B4', 'SR_B3', 'SR_B2']).multiply(0.0000275).add(-0.2)//.updateMask(mask)
print(TrueColor)
Export.image.toDrive({
  scale : 30,
  region: export_region,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  image: TrueColor,
  folder: 'global_et_example',
  description: station_name + 'True',
  fileNamePrefix: station_name + '_TrueColor_LC08_133033_20170716',
})
