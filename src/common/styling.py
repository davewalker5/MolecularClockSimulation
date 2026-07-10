from common.constants import DARK_THEME


def dark_theme_css() -> str:
    """Build app-specific CSS for the Field Notes dark theme.

    :return: CSS rules to inject into the Streamlit app.
    """
    colors = DARK_THEME
    # Streamlit's built-in theme handles the base palette; these rules tune widgets.
    return f"""
    <style>
    :root {{
      color-scheme: dark;
    }}

    .stApp {{
      background: {colors["page_bg"]};
      color: {colors["text"]};
    }}

    [data-testid="stSidebar"],
    [data-testid="stHeader"],
    [data-testid="stToolbar"] {{
      background: {colors["surface"]};
      border-color: {colors["border"]};
    }}

    [data-testid="stSidebar"] {{
      border-right: 1px solid {colors["border"]};
    }}

    h1, h2, h3, h4, h5, h6,
    .stMarkdown,
    .stText,
    label,
    p {{
      color: {colors["text"]};
    }}

    a {{
      color: {colors["link"]};
    }}

    a:hover {{
      color: {colors["link_hover"]};
    }}

    [data-testid="stMetric"],
    [data-testid="stCodeBlock"],
    [data-testid="stGraphVizChart"],
    div[data-baseweb="select"] > div,
    div[data-baseweb="input"] > div {{
      background: {colors["surface_elevated"]};
      border-color: {colors["border_strong"]};
      color: {colors["text"]};
    }}

    [data-testid="stMetric"] {{
      border: 1px solid {colors["border_strong"]};
      border-radius: 8px;
      padding: 0.75rem 0.9rem;
    }}

    [data-testid="stMetricLabel"],
    [data-testid="stCaptionContainer"],
    .st-emotion-cache-1wivap2 {{
      color: {colors["text_muted"]};
    }}

    .stTabs [data-baseweb="tab-list"] {{
      border-bottom: 1px solid {colors["border"]};
    }}

    .stTabs [data-baseweb="tab"] {{
      color: {colors["text_muted_strong"]};
    }}

    .stTabs [aria-selected="true"] {{
      color: {colors["white"]};
      border-bottom-color: {colors["link"]};
    }}

    .stButton > button,
    .stDownloadButton > button {{
      background: {colors["button"]};
      border-color: {colors["button"]};
      color: {colors["white"]};
      border-radius: 6px;
    }}

    .stButton > button:hover,
    .stDownloadButton > button:hover {{
      background: {colors["button_hover"]};
      border-color: {colors["button_hover"]};
      color: {colors["white"]};
    }}

    .stButton > button:disabled,
    .stDownloadButton > button:disabled {{
      background: {colors["row_stripe"]};
      border-color: {colors["border"]};
      color: {colors["text_subtle"]};
    }}
    </style>
    """
