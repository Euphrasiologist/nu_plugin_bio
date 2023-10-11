use nu_plugin::{serve_plugin, JsonSerializer};
use nu_plugin_bio::Bio;

fn main() {
    serve_plugin(&mut Bio {}, JsonSerializer {})
}
