/*****************************************************************
 * Function: fit
 * The function to fit the top/bottom edges using the given image. 
 * 
 * Input: image, type: ee.Image;
 * Input: X_name, band name for X axis, type: String/ee.String;
 * Input: Y_name, band name for Y axis, type: String/ee.String;
 * Input: X_from, the left of the range of variable X,  type: ee.Number, default: 0;
 * Input: X_to,   the right of the range of variable X, type: ee.Number, default: 1;
 * Input: X_step, the width of interval for variable X, type: ee.Number, default: 0.01;
 * 
 * Output: k and b for top/bottom edges, type: ee.Feature,
 * properties: [top_k, top_b, bottom_k, bottom_b];
*****************************************************************/
exports.fit = function(image, X_name, Y_name, X_from, X_to, X_step){
  X_from = (typeof X_from !== 'undefined') ?  X_from : 0;
  X_to   = (typeof X_to   !== 'undefined') ?  X_to   : 1;
  X_step = (typeof X_step !== 'undefined') ?  X_step : 0.01;
  
  var it = ee.List.sequence(X_from, X_to, X_step);
  
  var _top = it.map(function(n){
    var groups =  ee.Image(image).select(Y_name)
                    .updateMask(ee.Image(image).select(X_name).gte(ee.Number(n)))
                    .updateMask(ee.Image(image).select(X_name).lt(ee.Number(n).add(X_step)))
                    .reduceRegion({
                      geometry: ee.Image(image).geometry(),
                      reducer: ee.Reducer.max(),
                      maxPixels: 1e10
                    })
    var count = ee.Image(image).select(Y_name)
                  .updateMask(ee.Image(image).select(X_name).gte(ee.Number(n)))
                  .updateMask(ee.Image(image).select(X_name).lt(ee.Number(n).add(X_step)))
                  .reduceRegion({
                    geometry: ee.Image(image).geometry(),
                    reducer: ee.Reducer.count(),
                    maxPixels: 1e10
                  });
    return ee.Algorithms.If(
      ee.Number(ee.Dictionary(count).get(Y_name)).gt(5),
      ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(groups).get(Y_name)))
                      .set(X_name, ee.Number(n)),
      null);
  }, true);
  var tops = ee.List.sequence(X_from, X_to, ee.Number(X_step).multiply(5)).map(function(n){
    var groups = ee.FeatureCollection(ee.List(_top))
                   .filterMetadata(X_name, 'greater_than', ee.Number(n))
                   .filterMetadata(X_name, 'less_than', ee.Number(n).add(ee.Number(X_step).multiply(5)))
    return ee.Algorithms.If(
      groups.size().gte(3),
      ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [Y_name]})).get('median')))
                      .set(X_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [X_name]})).get('median'))),
      null);
  }, true);

  var _bottom = it.map(function(n){
    var groups =  ee.Image(image).select(Y_name)
                    .updateMask(ee.Image(image).select(X_name).gte(ee.Number(n)))
                    .updateMask(ee.Image(image).select(X_name).lt(ee.Number(n).add(X_step)))
                    .reduceRegion({
                      geometry: ee.Image(image).geometry(),
                      reducer: ee.Reducer.min(),
                      maxPixels: 1e10
                    })
    var count = ee.Image(image).select(Y_name)
                  .updateMask(ee.Image(image).select(X_name).gte(ee.Number(n)))
                  .updateMask(ee.Image(image).select(X_name).lt(ee.Number(n).add(X_step)))
                  .reduceRegion({
                    geometry: ee.Image(image).geometry(),
                    reducer: ee.Reducer.count(),
                    maxPixels: 1e10
                  });
    return ee.Algorithms.If(
      ee.Number(ee.Dictionary(count).get(Y_name)).gt(5),
      ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(groups).get(Y_name)))
                      .set(X_name, ee.Number(n)),
      null);
  }, true);
  var bottoms = ee.List.sequence(X_from, X_to, ee.Number(X_step).multiply(5)).map(function(n){
    var groups = ee.FeatureCollection(ee.List(_bottom))
                   .filterMetadata(X_name, 'greater_than', ee.Number(n))
                   .filterMetadata(X_name, 'less_than', ee.Number(n).add(ee.Number(X_step).multiply(5)))
    return ee.Algorithms.If(
      groups.size().gte(3),
      ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [Y_name]})).get('median')))
                      .set(X_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [X_name]})).get('median'))),
      null);
  }, true);
  
  
  var top   = ee.FeatureCollection(ee.List(tops)).reduceColumns({reducer : ee.Reducer.linearFit(), selectors : [X_name, Y_name]});
  var top_k = top.get("scale");
  var top_b = top.get("offset");
  
  var bottom   = ee.FeatureCollection(ee.List(bottoms)).reduceColumns({reducer : ee.Reducer.linearFit(), selectors : [X_name, Y_name]});
  var bottom_k = bottom.get("scale");
  var bottom_b = bottom.get("offset");
  
  return ee.Algorithms.If({
    condition: top_k,
    trueCase : ee.Feature(null).set('top_k', top_k)
                               .set('top_b', top_b)
                               .set('bottom_b', bottom_b)
                               .set('bottom_k', bottom_k),
    falseCase: ee.Feature(null)
  });
}

/*********************************************************************
 * Function: predict
 * The function to fit the dry/wet edges using the given image. 
 * 
 * Input: image,  type: ee.Image;
 * Input: para,   type: Strinee.Feature,
 *        properties: [top_k, top_b, bottom_k, bottom_b];
 * Input: X_name, band name for X axis, type: String/ee.String;
 * Input: Y_name, band name for Y axis, type: String/ee.String;
 * 
 * Output: remote sensing index from trapezoid model, type: ee.Image,
 * 1 means top, and 0 means bottom;
**********************************************************************/
exports.predict = function(image, para, X_name, Y_name){
  var top_k = ee.Number(ee.Feature(para).get('top_k'));
  var top_b = ee.Number(ee.Feature(para).get('top_b'));
  var bottom_k = ee.Number(ee.Feature(para).get('bottom_k'));
  var bottom_b = ee.Number(ee.Feature(para).get('bottom_b'));
  var X  = ee.Image(image).select(X_name);
  var Y  = ee.Image(image).select(Y_name);
  var rs = Y.subtract(X.multiply(bottom_k).add(bottom_b))
            .divide(X.multiply(top_k).add(top_b)
                     .subtract(X.multiply(bottom_k).add(bottom_b))
            );
  return rs;
}

/*****************************************************************
 * Function: fit_collection
 * The function to fit the top/bottom edges using the given image. 
 * 
 * Input: images, type: ee.ImageCollection;
 * Input: X_name, band name for X axis, type: String/ee.String;
 * Input: Y_name, band name for Y axis, type: String/ee.String;
 * Input: X_from, the left of the range of variable X,  type: ee.Number, default: 0;
 * Input: X_to,   the right of the range of variable X, type: ee.Number, default: 1;
 * Input: X_step, the width of interval for variable X, type: ee.Number, default: 0.01;
 * 
 * Output: k and b for top/bottom edges, type: ee.Feature,
 * properties: [top_k, top_b, bottom_k, bottom_b];
*****************************************************************/
exports.fit_collection = function(images, X_name, Y_name, X_from, X_to, X_step){
  X_from = (typeof X_from !== 'undefined') ?  X_from : 0;
  X_to   = (typeof X_to   !== 'undefined') ?  X_to   : 1;
  X_step = (typeof X_step !== 'undefined') ?  X_step : 0.01;
  
  var it = ee.List.sequence(X_from, X_to, X_step);
  
  var _top = it.iterate(function(n, feas){
    var groups  = ee.ImageCollection(images)
                    .map(function(image){
                      var _g =  ee.Image(image).select(Y_name)
                                  .updateMask(ee.Image(image).select(X_name).gte(ee.Number(n)))
                                  .updateMask(ee.Image(image).select(X_name).lt(ee.Number(n).add(X_step)))
                                  .reduceRegion({
                                    geometry: ee.Image(image).geometry(),
                                    reducer: ee.Reducer.max(),
                                    maxPixels: 1e10
                                  })
                      var count = ee.Image(image).select(Y_name)
                                    .reduceRegion({
                                      geometry: ee.Image(image).geometry(),
                                      reducer: ee.Reducer.count(),
                                      maxPixels: 1e10
                                    });
                      return ee.Algorithms.If(
                        ee.Number(ee.Dictionary(count).get(Y_name)).gt(5),
                        ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(_g).get(Y_name)))
                                        .set(X_name, ee.Number(n)),
                        null); 
                    }, true)
    return ee.FeatureCollection(feas).merge(groups);
  }, ee.FeatureCollection([]));
  var tops = ee.List.sequence(X_from, X_to, ee.Number(X_step).multiply(5)).map(function(n){
    var groups = ee.FeatureCollection(_top)
                   .filterMetadata(X_name, 'greater_than', ee.Number(n))
                   .filterMetadata(X_name, 'less_than', ee.Number(n).add(ee.Number(X_step).multiply(5)))
    return ee.Algorithms.If(
      groups.size().gte(3),
      ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [Y_name]})).get('median')))
                      .set(X_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [X_name]})).get('median'))),
      null);
  }, true);

  var _bottom = it.iterate(function(n, feas){
    var groups  = ee.ImageCollection(images)
                    .map(function(image){
                      var _g =  ee.Image(image).select(Y_name)
                                  .updateMask(ee.Image(image).select(X_name).gte(ee.Number(n)))
                                  .updateMask(ee.Image(image).select(X_name).lt(ee.Number(n).add(X_step)))
                                  .reduceRegion({
                                    geometry: ee.Image(image).geometry(),
                                    reducer: ee.Reducer.min(),
                                    maxPixels: 1e10
                                  })
                      var count = ee.Image(image).select(Y_name)
                                    .reduceRegion({
                                      geometry: ee.Image(image).geometry(),
                                      reducer: ee.Reducer.count(),
                                      maxPixels: 1e10
                                    });
                      return ee.Algorithms.If(
                        ee.Number(ee.Dictionary(count).get(Y_name)).gt(5),
                        ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(_g).get(Y_name)))
                                        .set(X_name, ee.Number(n)),
                        null); 
                    }, true)
    return ee.FeatureCollection(feas).merge(groups);
  }, ee.FeatureCollection([]));
  var bottoms = ee.List.sequence(X_from, X_to, ee.Number(X_step).multiply(5)).map(function(n){
    var groups = ee.FeatureCollection(_bottom)
                   .filterMetadata(X_name, 'greater_than', ee.Number(n))
                   .filterMetadata(X_name, 'less_than', ee.Number(n).add(ee.Number(X_step).multiply(5)))
    return ee.Algorithms.If(
      groups.size().gte(3),
      ee.Feature(null).set(Y_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [Y_name]})).get('median')))
                      .set(X_name, ee.Number(ee.Dictionary(groups.reduceColumns({reducer : ee.Reducer.median(), selectors : [X_name]})).get('median'))),
      null);
  }, true);
  
  
  var top   = ee.FeatureCollection(ee.List(tops)).reduceColumns({reducer : ee.Reducer.linearFit(), selectors : [X_name, Y_name]});
  var top_k = top.get("scale");
  var top_b = top.get("offset");
  
  var bottom   = ee.FeatureCollection(ee.List(bottoms)).reduceColumns({reducer : ee.Reducer.linearFit(), selectors : [X_name, Y_name]});
  var bottom_k = bottom.get("scale");
  var bottom_b = bottom.get("offset");
  
  return ee.Algorithms.If({
    condition: top_k,
    trueCase : ee.Feature(null).set('top_k', top_k)
                               .set('top_b', top_b)
                               .set('bottom_b', bottom_b)
                               .set('bottom_k', bottom_k),
    falseCase: ee.Feature(null)
  });
}
