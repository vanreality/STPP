//
// This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
//

import org.everit.json.schema.Schema
import org.everit.json.schema.loader.SchemaLoader
import org.everit.json.schema.ValidationException
import org.json.JSONObject
import org.json.JSONTokener
import org.json.JSONArray
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

class Schema {

    //
    // Resolve Schema path relative to main workflow directory
    //
    public static String getSchemaPath(workflow, schema_filename='schema.json') {
        return "${workflow.projectDir}/${schema_filename}"
    }

    
}