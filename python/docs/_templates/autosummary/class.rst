{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% if attributes %}
      .. autosummary::
         :toctree:
      {% for item in all_attributes %}
         {{ name }}.{{ item }}
      {% endfor %}
   {% endif %}

   {% if methods %}
      .. autosummary::
         :toctree:
      {% for item in methods %}
         {{ name }}.{{ item }}
      {% endfor %}
   {% endif %}