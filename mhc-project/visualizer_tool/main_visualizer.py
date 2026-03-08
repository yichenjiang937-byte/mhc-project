import os
import sys
import threading
import webbrowser


def get_script_dir() -> str:
    if getattr(sys, 'frozen', False):
        return getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    os.environ.setdefault('STREAMLIT_BROWSER_GATHER_USAGE_STATS', 'false')
    os.environ.setdefault('STREAMLIT_SERVER_HEADLESS', 'true')
    os.environ.setdefault('STREAMLIT_SERVER_PORT', '8501')
    os.environ.setdefault('STREAMLIT_SERVER_ADDRESS', '127.0.0.1')

    script_path = os.path.join(get_script_dir(), 'app_visualization.py')
    if not os.path.exists(script_path):
        raise FileNotFoundError(f'找不到 app_visualization.py: {script_path}')

    threading.Timer(1.5, lambda: webbrowser.open('http://127.0.0.1:8501')).start()

    from streamlit.web import bootstrap
    bootstrap.run(script_path, False, [], {})


if __name__ == '__main__':
    main()
