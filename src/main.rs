use nu_plugin::{serve_plugin, MsgPackSerializer};
use nu_plugin_bio::Bio;

fn main() {
    serve_plugin(&mut Bio {}, MsgPackSerializer {})
}