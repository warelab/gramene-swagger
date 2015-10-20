'use strict';

var SwaggerExpress = require('swagger-express-mw');
var express = require('express');
var app = express();
//var fs = require('fs');
//var yaml = require('js-yaml');

module.exports = app;

var config = {
  appRoot: __dirname // required config
};

// load the YAML-encoded schema. We will stream this out as JSON for swagger-ui
//var schema = yaml.load(fs.readFileSync(__dirname + '/api/swagger/swagger.yaml', 'utf8'));

SwaggerExpress.create(config, function (err, swaggerExpress) {
  if (err) { throw err; }

  // install middleware
  swaggerExpress.register(app);

  var port = process.env.PORT || 10010;
  app.listen(port);

  console.log('Listening on', port);

  // define routes
  app.get('/', function (req, res, next) { // redirect to /docs with the correct schema
    res.redirect('/docs?url=/swagger');
  });

  app.use('/docs', express.static('node_modules/swagger-ui/dist'));

  //app.get('/gramene.json', function (req, res, next) { // return top level info
  //  res.json(schema);
  //});

});

