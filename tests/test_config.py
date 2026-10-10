import pytest
from fakes import DT

from famescan.config import Config, FameClass, load_config, scans, to_toml


def test_default_durations_in_scans():
    # at 0.679 s per scan
    cfg = Config()
    assert scans(cfg.smoothing_window_s, DT) == scans(cfg.min_separation_s, DT) == 5
    assert cfg.min_width_s / DT == pytest.approx(1.3, abs=0.005)


def test_load_config(tmp_path):
    path = tmp_path / "setup.toml"
    path.write_text(
        'fame_table = "fames.csv"\n'
        'standard = "C17:0"\n'
        "noise_factor = 4\n"
        'report_fames = ["C16:0", "C18:1"]\n'
        "\n"
        "[classes.1]\n"
        "ecl_window = [0, 1.5]\n"
        "\n"
        "[classes.3]\n"
        "main_ions = { 79 = 40, 91 = 30 }\n"
        "ecl_window = [0, 2]\n"
    )
    cfg = load_config(path, standard_in_pct=False, standard=None)
    assert cfg.fame_table == str(tmp_path / "fames.csv")
    assert cfg.standard == "C17:0"
    assert cfg.noise_factor == 4.0 and isinstance(cfg.noise_factor, float)
    assert cfg.report_fames == ("C16:0", "C18:1")
    default = Config().classes
    assert cfg.classes[0] == default[0]
    assert cfg.classes[1] == FameClass(default[1].main_ions, (0.0, 1.5))  # main ions kept
    assert cfg.classes[3] == FameClass({79: 40.0, 91: 30.0}, (0.0, 2.0))
    assert cfg.standard_in_pct is False
    assert load_config(path, standard="").standard == ""


@pytest.mark.parametrize(
    ("toml", "error"),
    [
        ("smoothing_window = 5", "unknown setting.*smoothing_window"),
        ('report_fames = "C16:0"', "report_fames must be a list"),
        ("[classes.1]\nions = { 79 = 40 }", "classes.1: unknown setting"),
        ("[classes.3]\nmain_ions = { 79 = 40 }", "classes.3 needs main_ions and ecl_window"),
    ],
)
def test_load_config_errors(tmp_path, toml, error):
    path = tmp_path / "setup.toml"
    path.write_text(toml + "\n")
    with pytest.raises(ValueError, match=error):
        load_config(path)


def test_toml_roundtrip(tmp_path):
    cfg = Config(report_fames=("C16:0",), standard="", standard_in_pct=False)
    path = tmp_path / "setup.toml"
    path.write_text(to_toml(cfg))
    assert load_config(path) == cfg
    path.write_text(to_toml(Config()))
    assert load_config(path) == Config()
